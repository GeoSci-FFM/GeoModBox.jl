using Plots, Interpolations, Statistics
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.InitialCondition, GeoModBox.Tracers.TwoD
using Base.Threads
using Printf

@doc raw"""
    Advection_2D_ResTest()

...

"""
function Advection_2D_ResTest()

@printf("Running on %d thread(s)\n", nthreads())

nrnxny      =   10
Scheme      =   ["upwind","slf","semilag","tracers"]
ns          =   size(Scheme,1)
@show ns

# Statistical Parameter ============================================== #
St      = (
    Δ           =   zeros(size(Scheme,1),nrnxny),    
    nxny        =   zeros(size(Scheme,1),nrnxny),
    Tmax        =   zeros(size(Scheme,1),nrnxny),
    Tmean       =   zeros(size(Scheme,1),nrnxny),    
)
# -------------------------------------------------------------------- #
# Define Initial Condition =========================================== #
# Temperature - 
#   1) circle, 2) gaussian, 3) block
# Velocity - 
#   1) RigidBody, 2) ShearCell
Ini         =   (T=:circle,V=:RigidBody,) 
# -------------------------------------------------------------------- #
for m = 1:ns # Loop over advection schemes
    # Define Numerical Scheme ======================================== #
    FD          =   (Method     = (Adv=Scheme[m],),)    
    @printf("Advection Scheme: %s\n ",string(FD.Method.Adv))
    # ---------------------------------------------------------------- #
    # Plot constants ================================================= #
    Pl  =   (
        inc         =   5,
        sc          =   1.0e9,
        Minc        =   1, 
        Msz         =   0.2,
    )
    # ---------------------------------------------------------------- #
    # Model Constants ================================================ #
    M   =   (
        xmin    =   0.0,
        xmax    =   1.0,
        ymin    =   0.0,
        ymax    =   1.0,
    )
    # ---------------------------------------------------------------- #
    for l = 1:nrnxny # Loop over differnet resolutions
        # Numerical Constants ======================================== #
        NC  =   (
            x       =   l*20,       # Number of horizontal centroids
            y       =   l*20,       # Number of vertical centroids
        )
        display(string("nx = ",NC.x,", ny = ",NC.y))
        NV =   (
            x       =   NC.x + 1,   # Number of horizontal vertices
            y       =   NC.y + 1,   # Number of vertical vertices
        )
        Δ   =   (
            x   =   (abs(M.xmin)+M.xmax)/NC.x,
            y   =   (abs(M.ymin)+M.ymax)/NC.y,
        )
        # ------------------------------------------------------------ #
        # Grid ======================================================= #
        x   =   (
            c       =   LinRange(M.xmin + Δ.x/2.0, M.xmax - Δ.x/2.0, NC.x),
            ce      =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
            v       =   LinRange(M.xmin, M.xmax , NV.x)    
        )
        y       = (
            c       =   LinRange(M.ymin + Δ.y/2.0, M.ymax - Δ.y/2.0, NC.y),
            ce      =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
            v       =   LinRange(M.ymin, M.ymax, NV.y),    
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
        # ------------------------------------------------------------ #
        # Animationsettings =0======================================== #
        path        =   string("./examples/AdvectionEquation/Results/")
        anim        =   Plots.Animation(path, String[] )
        filename    =   string("2D_advection_",Ini.T,"_",Ini.V,
                                "_",FD.Method.Adv,"_",NC.x,"_",NC.y,
                                "_nth_",nthreads())
        save_fig    =   -1
        # ------------------------------------------------------------ #
        # Array Initialization ======================================= #
        D       =   (
            T       =   zeros(Float64,(NC.x,NC.y)),
            T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
            T_exo   =   zeros(Float64,(NC.x+2,NC.y+2)),
            vx      =   zeros(Float64,(NV.x,NV.y+1)),
            vy      =   zeros(Float64,(NV.x+1,NV.y)),    
            vxc     =   zeros(Float64,(NC.x,NC.y)),
            vyc     =   zeros(Float64,(NC.x,NC.y)),
            vc      =   zeros(Float64,(NC.x,NC.y)),
            wt      =   zeros(Float64,(NC.x,NC.y)),
            Tmax    =   [0.0],
            Tmin    =   [0.0],
            Tmean   =   [0.0],
        )
        # Initial Condition ========================================== #
        # Temperature ---
        IniTemperature!(Ini.T,M,NC,Δ,D,x,y)
        if FD.Method.Adv == "slf"
            D.T_exo    .=  D.T_ex
        end
        # Velocity ---
        IniVelocity!(Ini.V,D,NV,Δ,M,x,y)        
        # Get the velocity on the centroids ---
        @threads for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # ------------------------------------------------------------ #
        # Time ======================================================= #
        T   =   ( 
            tmax    =   [0.0],  
            Δfac    =   1.0,    # Courant time factor, i.e. dtfac*dt_courant
            Δ       =   [0.0],
        )
        T.tmax[1]   =   π*((M.xmax-M.xmin)-Δ.x)/maximum(D.vc)   # t = U/v [ s ]
        T.Δ[1]      =   T.Δfac * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        nt          =   ceil(Int,T.tmax[1]/T.Δ[1])
        # ------------------------------------------------------------ #
        # Tracer Advection =========================================== #
        if FD.Method.Adv == "tracers"
            # Tracer Initialization ---
            nmx,nmy     =   3,3
            noise       =   1
            nmark       =   nmx*nmy*NC.x*NC.y
            Aparam      =   :thermal
            MPC         =   (
                c               =   zeros(Float64,(NC.x,NC.y)),
                th              =   zeros(Float64,(nthreads(),NC.x,NC.y)),
                min             =   zeros(Float64,nt),
                max             =   zeros(Float64,nt),
                mean            =   zeros(Float64,nt),
            )
            MPC1        = (
                PG_th   =   [similar(D.T) for _ = 1:nthreads()], # per thread
                wt_th   =   [similar(D.wt) for _ = 1:nthreads()], # per thread
            )
            MPC     =   merge(MPC,MPC1)
            Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise)
            # RK4 weights ---
            rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
            rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
            # Interpolate on centroids ---
            @threads for k = 1:nmark
                Ma.T[k] =   FromCtoM(D.T_ex, k, Ma, x, y, Δ, NC)
            end
            # Count marker per cell ---
            CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,1)
        end
        # ------------------------------------------------------------ #
        # Visualize initial condition -------------------------------- #
        if FD.Method.Adv == "tracers"
            #p = scatter(Ma.x[1:Pl.Minc:end],Ma.y[1:Pl.Minc:end], 
            #        zcolor=Ma.T[1:Pl.Minc:end],color=:thermal,
            #        markersize=Pl.Msz,legend=false,colorbar=true, 
            #        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
            #        layout=(1,2),subplot=1)
            p = heatmap(x.c,y.c,(D.T./D.Tmax)',color=:thermal, 
                    aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
                    ylims=(M.ymin, M.ymax),clims=(0.5, 1.0),
                    colorbar=true,layout=(1,2),subplot=1)
            quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
                    y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                    quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                            D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                    color="white",layout=(1,2),subplot=1)
            heatmap!(p,x.c,y.c,MPC.c',color=:inferno, 
                    aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
                    colorbar=true,clims=(0.0, 18.0),title=:"Marker per cell",
                    layout=(1,2),subplot=2)
        else
            p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
                    color=:thermal, colorbar=true, aspect_ratio=:equal, 
                    xlabel="x", ylabel="z", 
                    title="Temperature", 
                    xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
                    clims=(0.5, 1.0))
            quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                    quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                            D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                    color="white")
        end
        if save_fig == 1
            Plots.frame(anim)
        elseif save_fig == 0
            display(p)
        end
        # ------------------------------------------------------------ #
        # Time Loop ================================================== #
        for i=2:nt
            #@printf("Time step: #%04d\n ",i) 

            if FD.Method.Adv == "upwind"
                upwindc2D!(D.T,D.T_ex,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
            elseif FD.Method.Adv == "slf"
                slfc2D!(D.T,D.T_ex,D.T_exo,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
            elseif FD.Method.Adv == "semilag"
                semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,[],[],x,y,T.Δ[1])
            elseif FD.Method.Adv == "tracers"
                # Advect tracers ---
                AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
                CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,i)
                
                # Interpolate temperature from tracers to grid ---
                Markers2Cells(Ma,nmark,MPC.PG_th,D.T,MPC.wt_th,D.wt,x,y,Δ,Aparam)        
            end
            
            display(string("ΔT = ",((maximum(filter(!isnan,D.T))-D.Tmax[1])/D.Tmax[1])*100))

            # Plot Solution ---
            if mod(i,10) == 0 || i == nt
                if FD.Method.Adv == "tracers"
                    #p = scatter(Ma.x[1:Pl.Minc:end],Ma.y[1:Pl.Minc:end], 
                    #        zcolor=Ma.T[1:Pl.Minc:end],color=:thermal,
                    #        markersize=Pl.Msz,legend=false,colorbar=true, 
                    #        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
                    #        layout=(1,2),subplot=1)
                    p = heatmap(x.c,y.c,(D.T./D.Tmax)',color=:thermal, 
                            aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
                            ylims=(M.ymin, M.ymax),clims=(0.5, 1.0),
                            colorbar=true,layout=(1,2),subplot=1)
                    quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
                            y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                            quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                                    D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                            color="white",layout=(1,2),subplot=1)
                    heatmap!(p,x.c,y.c,MPC.c',color=:inferno, 
                            aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
                            colorbar=true,clims=(0.0, 18.0),title=:"Marker per cell",
                            layout=(1,2),subplot=2)
                else
                    p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
                            color=:thermal, colorbar=true, aspect_ratio=:equal, 
                            xlabel="x", ylabel="z", 
                            title="Temperature", 
                            xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
                            clims=(0.5, 1.0))
                    quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
                                y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                                quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                                        D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                            color="white")
                end
                if save_fig == 1
                    Plots.frame(anim)
                elseif save_fig == 0
                    display(p)                        
                end
            end
                
        end # End Time loop    
        # Save Animation ============================================= #
        if save_fig == 1
            # Write the frames to a GIF file
            Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
            foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
        elseif save_fig == 0
            display(plot(p))
        end
        # Statistical Values for Each Scheme and Resolution ---
        St.Δ[m,l]       =   abs((maximum(filter(!isnan,D.T))-D.Tmax[1])/D.Tmax[1])*100
        St.nxny[m,l]    =   1/NC.x/NC.y
        St.Tmax[m,l]    =   maximum(filter(!isnan,D.T))
        St.Tmean[m,l]   =   mean(abs.(filter(!isnan,D.T)))
        # ------------------------------------------------------------ #

    end # End resolution loop

end # End method loop

q   =   plot(0,0,layout=(1,3))
for m=1:ns    
    plot!(q,St.nxny[m,:],St.Δ[m,:],
                marker=:circle,markersize=3,label=Scheme[m],
                xaxis=:log,yaxis=:log,
                xlims=(minimum(St.nxny), maximum(St.nxny)), 
                ylims=(1e-15, 1e3), 
                xlabel="1/nx/ny",ylabel="ΔT[%]",layout=(1,3),
                subplot=1)
    plot!(q,St.nxny[m,:],St.Tmax[m,:],
                marker=:circle,markersize=3,label="",
                xaxis=:log,yaxis=:log,
                xlims=(minimum(St.nxny), maximum(St.nxny)), 
                ylims=(1e2, 1e5),
                xlabel="1/nx/ny",ylabel="T_{max}",
                subplot=2)
        #plot!(q,1/maximum(St.nxny[1,:]):1e-4:1/minimum(St.nxny[1,:]),
        #        St.Tanamax .* ones(1,length(1/maximum(St.nxny[1,:]):1e-4:1/minimum(St.nxny[1,:]))),
        #        linecolor=:black,linestyle=:dash)
    plot!(q,St.nxny[m,:],St.Tmean[m,:],
                marker=:circle,markersize=3,label="",
                xaxis=:log,yaxis=:log,
                xlims=(minimum(St.nxny), maximum(St.nxny)), 
                ylims=(1e2, 1e4), 
                xlabel="1/nx/ny",ylabel="⟨T⟩",
                subplot=3)
    display(q)
end
# --------------------------------------------------------------------- #
# Save Final Figure =================================================== #
savefig(q,string("./examples/AdvectionEquation/",
                    "Results/2D_advection_",Ini.T,"_",
                    Ini.V,"_ResTest.png"))
# --------------------------------------------------------------------- #

end # Function end

Advection_2D_ResTest()