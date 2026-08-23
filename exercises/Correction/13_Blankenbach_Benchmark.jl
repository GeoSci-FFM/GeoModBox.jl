using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Scaling
using Statistics, LinearAlgebra
using Printf, LaTeXStrings, Measures

function BlankenbachBenchmark()
# Define Initial Condition ============================================== #
# Temperature - 
#   1) circle, 2) gaussian, 3) block, 4) linear, 5) lineara
# !!! Gaussian is not working!!! 
Ini         =   (T=:blankenbach,)
# ----------------------------------------------------------------------- #
# Plot Settings ========================================================= #
Pl  =   (
    qinc        =   10,
    qsc         =   1.0e-3 # 1.0e-3, 4.0e-4, 5.0e-5
)
# Animation Settings ==================================================== #
path        =   string("./exercises/Correction/Results/")
anim        =   Plots.Animation(path, String[] )
save_fig    =   1
# ----------------------------------------------------------------------- #
# Benchmark Values ====================================================== # 
# Taken from Gerya (2019), Introduction to numerical geodynamic modelling
B       =   (
    Nr      =   1,
    # Nusselt Number at the top
    Nu      =   [4.8844,10.534,21.972,10.066,6.9229],   
    # Root Mean Square Velocity-scaled
    Vrms    =   [42.865,193.21,833.99,480.43,171.76],   
    # Non-dimensional temperature gradients in the model corner
    q1      =   [8.0593,19.079,45.964,17.531,18.484],
    q2      =   [0.5888,0.7228,0.8772,1.0085,0.1774],
    q3      =   [8.0593,19.079,45.964,26.809,14.168], 
    q4      =   [0.5888,0.7228,0.8772,0.4974,0.6177],
    # Local minimum along the central vertical temperature profile 
    Tmin    =   [0.4222,0.4284,0.4322,0.7405,0.3970],
    ymin    =   [0.2249,0.1118,0.0577,0.0623,0.1906],
    # Local maximum algon the central vertical temperature profile
    Tmax    =   [0.5778,0.5716,0.5678,0.8323,0.5758],
    ymax    =   [0.7751,0.8882,0.9423,0.8243,0.7837],
)
# ----------------------------------------------------------------------- #
# Geometry ============================================================== #
M   =   Geometry(
    xmin    =   0.0,                #   [ m ] 
    xmax    =   1000e3,             #   [ m ]
    ymin    =   -1000e3,            #   [ m ]
    ymax    =   0.0,                #   [ m ]
)
# ----------------------------------------------------------------------- #
# Referenzparameter ===================================================== #
P   =   Physics(
    g       =   10.0,               #   Graviational Acceleration [m/s^2]
    ρ₀      =   4000.0,             #   Reference Density [kg/m^3]
    k       =   5.0,                #   Thermal Conductivity [ W/m/K ]
    cp      =   1250.0,             #   Heat capacity [ J/kg/K ]
    α       =   2.5e-5,             #   THermal Expansion [ K^-1 ]
    η₀      =   1e23,               #   Reference Viscosity [ Pa*s ]
    ΔT      =   1000.0,             #   Temperature Difference
    # If Ra < 0, Ra will be calculated from the reference parameters.
    # If Ra is defined, the reference viscosity will be adjusted such that
    # the scaling parameters result in the given Ra.
    Ra      =   -9999,              #   Rayleigh number
)
# ----------------------------------------------------------------------- #
# Define Scaling Constants ============================================== # 
S   =   ScalingConstants!(M,P)
# ----------------------------------------------------------------------- #
# Numerical Grid ======================================================== #
NC  =   (
    x   =   50,
    y   =   50,
)
NV      =   (
    x   =   NC.x + 1,
    y   =   NC.y + 1,
)
Δ       =   GridSpacing(
    x   =   (M.xmax - M.xmin)/NC.x,
    y   =   (M.ymax - M.ymin)/NC.y,
)
# ----------------------------------------------------------------------- #
# Data Arrays =========================================================== #
D       =   DataFields(
    Q       =   zeros(Float64,(NC...)),
    T       =   zeros(Float64,(NC...)),
    T0      =   zeros(Float64,(NC...)),
    T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
    T_ex0   =   zeros(Float64,(NC.x+2,NC.y+2)),
    Told_ex =   zeros(Float64,(NC.x+2,NC.y+2)),
    ρ       =   ones(Float64,(NC...)),
    vx      =   zeros(Float64,(NV.x,NV.y+1)),
    vy      =   zeros(Float64,(NV.x+1,NV.y)),    
    Pt      =   zeros(Float64,(NC...)),
    vxc     =   zeros(Float64,(NC...)),
    vyc     =   zeros(Float64,(NC...)),
    vxco    =   zeros(Float64,(NC...)),
    vyco    =   zeros(Float64,(NC...)),
    vc      =   zeros(Float64,(NC...)),
    ΔTtop   =   zeros(Float64,NC.x),
    ΔTbot   =   zeros(Float64,NC.x),
)
# ----------------------------------------------------------------------- #
# Needed for the defect correction solution ---
divV        =   zeros(Float64,NC...)
ε           =   (
    xx      =   zeros(Float64,NC...), 
    yy      =   zeros(Float64,NC...), 
    xy      =   zeros(Float64,NV...),
)
τ           =   (
    xx      =   zeros(Float64,NC...), 
    yy      =   zeros(Float64,NC...), 
    xy      =   zeros(Float64,NV...),
)
# Residuals ---
Fm     =    (
    x       =   zeros(Float64,NV.x, NC.y), 
    y       =   zeros(Float64,NC.x, NV.y)
)
FPt         =   zeros(Float64,NC...)
# Heat Flux calculations ---
Tv1         =   zeros(NV.x,1)
Tv2         =   zeros(NV.x,1)
dTdy        =   zeros(NV.x,1)
# ----------------------------------------------------------------------- #
# Time parameters ======================================================= #
T   =   TimeParameter(
    tmax    =   1000000.0,          #   [ Ma ]
    Δfacc   =   0.9,                #   Courant time factor
    Δfacd   =   0.9,                #   Diffusion time factor
    itmax   =   30000,              #   Maximum iterations
)
T.tmax      =   T.tmax*1e6*T.year   #   [ s ]
T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
T.Δd        =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))

T.Δ         =   minimum([T.Δd,T.Δc])
# Statistics ------------------------------------------------------------ #
Time            =   zeros(T.itmax)
Nu              =   zeros(T.itmax)
meanV           =   zeros(T.itmax)
meanT           =   zeros(T.itmax)
epsV_history    =   fill(NaN, T.itmax)
find            =   0
count           =   0
epsV0           =   0
final_step      =   0 
# ----------------------------------------------------------------------- #
# Scaling laws ========================================================== #
ScaleParameters!(S,M,Δ,T,P,D)
# ----------------------------------------------------------------------- #
# Coordinates =========================================================== #
x       =   (
    c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
    ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
    v   =   LinRange(M.xmin,M.xmax,NV.x),
)
y       =   (
    c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
    ce  =   LinRange(M.ymin - Δ.y/2.0, M.ymax + Δ.y/2.0, NC.y+2),
    v   =   LinRange(M.ymin,M.ymax,NV.y),
)
x1      =   (
    c2d     =   x.c .+ 0*y.c',
    v2d     =   x.v .+ 0*y.v', 
    vx2d    =   x.v .+ 0*y.ce',
    vy2d    =   x.ce .+ 0*y.v',
)
x   =   merge(x,x1)
y1      =   (
    c2d     =   0*x.c .+ y.c',
    v2d     =   0*x.v .+ y.v',
    vx2d    =   0*x.v .+ y.ce',
    vy2d    =   0*x.ce .+ y.v',
)
y   =   merge(y,y1)
# ----------------------------------------------------------------------- #
# Initial Condition ===================================================== #
# Temperature ------
IniTemperature!(Ini.T,M,NC,D,x,y;Tb=P.Tbot,Ta=P.Ttop)
# ----------------------------------------------------------------------- #
# Boundary Conditions =================================================== #
# Temperature ------
TBC     = (
    type    = (W=:Neumann, E=:Neumann,N=:Dirichlet,S=:Dirichlet),
    val     = (W=zeros(NC.y),E=zeros(NC.y),
                    N=P.Ttop.*ones(NC.x),S=P.Tbot.*ones(NC.x)))
# Velocity ------
VBC     =   (
    type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
    val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vyS=0.0,vyN=0.0,vxW=0.0,vxE=0.0),
)
# ----------------------------------------------------------------------- #
# Rayleigh Number ======================================================= #
P.Ra        =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.η₀/P.κ

filename    =   string("13_Blankenbach_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,
                        "_",Ini.T)
# ----------------------------------------------------------------------- #
# Linear System of Equations ============================================ #
# Momentum Conservation Equation (MCE) ------
niterM  =   50
atolM   =   1e-8        #   Absolute tolerance
rtolM   =   1e-5        #   Relative tolerance; r = atolM0/atolM
RM      =   0.0         #   Initialize absolute residual    
RMrel   =   0.0         #   Initialize relative residual 
off     =   [   NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y ] # Pt
Num     =    (
    Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
    Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
    Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
)
ndofM   =   maximum(Num.Pt)
KM      =   ExtendableSparseMatrix(ndofM,ndofM)      
δx      =   zeros(ndofM)
F       =   zeros(ndofM)
# Assemble Matrix for momentum equation ---
# Optimization outside time loop since viscosity is constant ---
KM      =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
KMfac   =   lu(KM.cscmatrix)
# Energy Conservation Equation (ECE) -------------------------------- #
niterT  =   10
ϵT      =   1e-12
RE      =   0.0
ndofT   =   maximum(Num.T)     
KT      =   ExtendableSparseMatrix(ndofT,ndofT)
δT      =   zeros(maximum(Num.T))
RT      =   zeros(Float64,NC...)
∂2T     =   (∂x2=zeros(NC.x, NC.y), ∂y2=zeros(NC.x, NC.y),
                ∂x20=zeros(NC.x, NC.y), ∂y20=zeros(NC.x, NC.y))
# ----------------------------------------------------------------------- #
# Time Loop ============================================================= #
for it = 1:T.itmax
    R0      =   0.0
    # Reduce screen output ---
    verbose_step    =   mod(it, 400) == 0 || it == 1 || final_step == 1
    if it>1
        Time[it]  =   Time[it-1] + T.Δ
    end
    verbose_step && @printf(
        "Time step: #%04d, Time [non-dim]: %04e\n",
        it, Time[it]
    )
    # ------ MCE ------
    verbose_step && @printf("---Momentum Calculation ---\n")
    if it == 1
        D.vx[:,2:end-1]     .=  0.0
        D.vy[2:end-1,:]     .=  0.0
        D.Pt                .=  0.0
    end
    # Residual Calculation ------
    @. D.ρ  =   -P.Ra*D.T
    for iter = 1:niterM
        Residuals2Dc!(D,VBC,ε,τ,divV,Δ,1.0,1.0,Fm,FPt)
        F[Num.Vx]   .=  Fm.x
        F[Num.Vy]   .=  Fm.y
        F[Num.Pt]   .=  FPt
        RM          =   norm(F)/sqrt(length(F))
        if iter == 1
            R0 = max(RM, eps())
        end
        RMrel       =   RM/R0
        if verbose_step
            @printf("   MCE %2d: ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n",iter,RM,RMrel)
        end
        (RM < atolM || RM/R0 < rtolM) && break
        # Solving Linear System of Equations ------
        δx      .=   -(KMfac \ F)
        # Update unknown variables ------
        D.vx[:,2:end-1]     .+=  δx[Num.Vx]
        D.vy[2:end-1,:]     .+=  δx[Num.Vy]
        D.Pt                .+=  δx[Num.Pt]
    end
    # ======
    # Calculate Velocity on the Centroids ------
    for i = 1:NC.x
        for j = 1:NC.y
            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
        end
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
    # Nusselt Number ==================================================== #
    # Grid structure at the surface ---
    #   o - Centroids
    #   x - Vertices 
    #   □ - Ghost Nodes
    #
    #   □          □           □            □
    #   
    #        x --------- x --------- x
    #        |           |           |
    #   □    |     o     |     o     |      □ 
    #        |           |           |
    #        x --------- x --------- x      
    #        |           |           |
    #   □    |     o     |     o     |      □
    # --- 
    # Get temperature at the vertices 
    @. Tv1  =   (D.T_ex[1:end-1,end] + D.T_ex[2:end,end] + 
                    D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1])/4
    @. Tv2  =   (D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1] + 
                    D.T_ex[1:end-1,end-2] + D.T_ex[2:end,end-2])/4
    # Calculate temperature gradient --- 
    @. dTdy =   -(Tv1 - Tv2)/Δ.y
    # Calculate Nusselt number ---
    Nu[it]      =   0.0
    # Trapezoidal integration -
    for i = 1:NV.x
        if i == 1 || i == NV.x
            afac = 1
        else
            afac = 2
        end
        Nu[it]     += afac * dTdy[i]
    end
    Nu[it]     *=   Δ.x/2
    # @show Nu[it]
    Nu[it]      /=  M.xmax - M.xmin
    # Mean Temperature ---
    meanT[it]   =   mean(D.T)
    # Root Mean Square Velocity ---
    meanV[it] = sqrt(mean(D.vxc.^2 .+ D.vyc.^2))
    # ------------------------------------------------------------------- #
    # Calculate time stepping =========================================== #
    T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
            (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd        =   T.Δfacd * (1.0 / (2.0 *(1.0/Δ.x^2 + 1/Δ.y^2)))
    T.Δ         =   minimum([T.Δd,T.Δc])
    if Time[it] + T.Δ >= T.tmax
        T.Δ        = T.tmax - Time[it]
        final_step = 1
    end
    # Plot ============================================================== #
    if verbose_step
        p = heatmap(x.c,y.c,D.T',
                xlabel= L"x",ylabel= L"y",colorbar=true,
                color=cgrad(:lajolla),
                aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
                ylims=(M.ymin, M.ymax), 
                size=(1100,900),dpi = 300,
                guidefontsize = 28, tickfontsize = 18,
                titlefontsize = 28, top_margin=5mm)
        contour!(p,x.c,y.c,D.T',lw=1,
                        color="white",cbar=false,
                        alpha=0.5)
        quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                        D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
                la=0.5,color="black",
                alpha=0.5)
        if save_fig == 1
            Plots.frame(anim)
        elseif save_fig == 0
            display(p)
        end
    end
    # ------------------------------------------------------------------- #
    # Advection ========================================================= #
    if it == 1
        @. D.vxco   =   D.vxc
        @. D.vyco   =   D.vyc
    end
    semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,D.vxco,D.vyco,x,y,T.Δ)
    @. D.vxco  =   D.vxc
    @. D.vyco  =   D.vyc
    # ------------------------------------------------------------------- #
    # Diffusion ========================================================= #
    verbose_step &&  @printf("---Energy Calculation---\n")
    # Update temperature field --- 
    @. D.T0 =   D.T
    # Assemble linear system ---
    KT      =   AssembleMatrix2Dc(1.0, TBC, Num, NC, Δ, T.Δ;C=0.5)
    # Solve for temperature correction: Cholesky factorisation
    KTc      =   cholesky(KT.cscmatrix)
    for iter = 1:niterT
        # Evaluate residual
        ComputeResiduals2Dc!( RT, D.T, D.T_ex, D.T0, D.T_ex0, ∂2T, 
                1.0, TBC, Δ, T.Δ; C = 0.5 )
        RE      =   norm(RT)/sqrt(length(RT))
        if verbose_step
                @printf("   ECE %2d: ||RE|| = %1.4e\n", iter, RE)
        end
        RE < ϵT ? break : nothing
        # Solve for temperature correction: Back substitutions
        δT .= -(KTc\RT[:])
        # Update temperature
        @. D.T += δT[Num.T]
    end
    # Update extended field for advection scheme --- 
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #
    # Check break ======================================================= #
    # If the maximum time is reached or if the models reaches steady 
    # state the time loop is stoped! 
    if Time[it] > 0.05
        epsC    =   0.001
        ind = searchsortedfirst(@view(Time[1:it]), Time[it] - 0.05)
        epsV    =   std(@view meanV[ind:it])
        if count == 0
            epsV0   = max(epsV, eps(Float64))
            count   += 1
        end
        epsV        /= epsV0
        epsV_history[it] = epsV
        find    =   it
        verbose_step && @printf("ε_Vr = %g, ε_Cr = %g \n",epsV,epsC)
        if Time[it] >= T.tmax
            @printf("Maximum time reached!\n")
            find    =   it
            break
        elseif (epsV <= epsC)
            @printf("Convection reaches steady state at %i!\n",it)
            find    =   it
            break
        end
    end
    # --------------------------------------------------------------- #
    verbose_step && @printf("\n")
end
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
    foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
end
valid   =   findall(isfinite, epsV_history)
    k   = plot(
            valid,
            log10.(epsV_history[valid]),
            xlabel = "it",
            ylabel = "log₁₀(εᵥ/εᵥ₀)",
            label = "",
            markershape = :circle,
            markercolor = :black,
    )
# Save final figure ===================================================== #
p2 = heatmap(x.c,y.c,D.T',
        xlabel= L"x",ylabel= L"y",colorbar=true,
        color=cgrad(:lajolla),
        aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
        ylims=(M.ymin, M.ymax), 
        size=(1100,900),dpi = 300,
        guidefontsize = 28, tickfontsize = 18,
        titlefontsize = 28, top_margin=5mm)
contour!(p2,x.c,y.c,D.T',lw=1,
                color="white",cbar=false,
                alpha=0.5)
quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
        quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
        la=0.5,color="black",
        alpha=0.5)
if save_fig == 1
    savefig(p2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Final_Stage_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,"_it_",find,"_",
            Ini.T,".png"))
    savefig(k,string("./exercises/Correction/Results/13_BlankenbachBenchmark_iterations_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,
            "_",Ini.T,".png"))
elseif save_fig == 0
    display(p2)
    display(k)
end
# ----------------------------------------------------------------------- #
# Plot time serieses ==================================================== #
q2  =   plot(layout=(2,1),size=(1200,900),dpi = 300)
plot!(q2,Time[1:find],Nu[1:find],
            xlabel= "", ylabel= L"Nu",label="",
            xformatter = _ -> "",
            ylims=(0,25),xlims=(0,Time[find]),
            guidefontsize = 14, tickfontsize = 14,
            titlefontsize = 14,grid = false,
                layout=(2,1),subplot=1)
plot!(q2,Time[1:find],B.Nu[B.Nr].*ones(find,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
                layout=(2,1),subplot=1)
plot!(q2,Time[1:find],meanV[1:find],
            xlabel= L"Time\ [\ non-dim\ ]", ylabel= L"V_{RMS}",label="",
            xformatter =:auto,
            guidefontsize = 14, tickfontsize = 14,
            titlefontsize = 14,grid = false,
            ylims=(0,8000),xlims=(0,Time[find]),
            layout=(2,1),subplot=2)
plot!(q2,Time[1:find],B.Vrms[B.Nr].*ones(find,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
            layout=(2,1),subplot=2)
if save_fig == 1
    savefig(q2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_TimeSeries_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
elseif save_fig == 0
    display(q2)
end
# ----------------------------------------------------------------------- #
# Vertical temperature profiles ========================================= #
ind1    =   Int64(floor((M.xmax-M.xmin)/2/Δ.x)+1)
ind2    =   Int64(floor((M.xmax-M.xmin)/2/Δ.x))
Tprof   =   (D.T[ind1,:]+D.T[ind2,:])/2
q3  =   plot(Tprof,y.c,
        xlabel= L"T",ylabel= L"y",label=false,
        size=(1100,900),dpi=300,grid = false,
        ylims=(M.ymin,M.ymax),xlims=(0.0,1.0),
        guidefontsize = 28, tickfontsize = 18,
        titlefontsize = 28,right_margin=5mm,
        left_margin=5mm)
scatter!(q3,(B.Tmin[B.Nr],-(1-B.ymin[B.Nr])),
            markershape=:rect,markersize=6,
            markercolor=:black,label=false)
scatter!(q3,(B.Tmax[B.Nr],-(1-B.ymax[B.Nr])),
            markershape=:rect,markersize=6,
            markercolor=:black,label=false)
if save_fig == 1
    savefig(q3,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Profiles_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
elseif save_fig == 0
    display(q3)
end
# ----------------------------------------------------------------------- #
end

start=time()
BlankenbachBenchmark()
stop=time()
@printf("Runtime: %g s",stop-start)