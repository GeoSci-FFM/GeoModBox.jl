using Plots
using ExtendableSparse
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD

function main()

    # Plot Settings ----------------------------------------------------- #
    Pl  =   (
        qinc    =   5,
        qsc     =   100*(60*60*24*365.25)*1e2
    )
    # ------------------------------------------------------------------- #
    # Geometry ---------------------------------------------------------- #
    M       =   (
        xmin    =   0.0,
        xmax    =   500.0e3,    # [ m ]
        ymin    =   -500.0e3,   # [ m ]
        ymax    =   0.0,
    )
    # -------------------------------------------------------------------- #
    # Grid --------------------------------------------------------------- #
    NC      =   (
        x   =   50, 
        y   =   50,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   (
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    x       =   (
        c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v   =   LinRange(M.xmin,M.xmax,NV.x),
    )
    y       =   (
        c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
        ce  =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
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
    # -------------------------------------------------------------------- #
    # Physics ------------------------------------------------------------ #
    g       =   9.81            #   Gravitational Acceleration 
    # Viscosity ---
    η₀      =   1.0e21          #   Background
    η₁      =   1.0e23          #   Anomaly
    η       =   [η₀,η₁]
    # Density ---
    ρ₀      =   3200.0          #   Background 
    ρ₁      =   3300.0          #   Anomaly
    ρ       =   [ρ₀,ρ₁]    
    
    phase   =   [0,1]
    # ------------------------------------------------------------------- #
    # Allocation -------------------------------------------------------- #
    D   =   (
        vx  =   zeros(NV.x,NC.y+2),
        vy  =   zeros(NC.x+2,NV.y),
        Pt  =   zeros(NC...),
        ρ   =   zeros(NC...),
        ηs  =   (c=zeros(NC...), v=zeros(NV...)),
        vxc =   zeros(NC...),
        vyc =   zeros(NC...),
        vc  =   zeros(NC...),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ------------------------------------------------- #
    type    =:block
    IniPhase!(type,D,M,x,y,NC;ρ,η)
    p   =   heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                xlabel="x[km]",ylabel="y[km]",
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                layout=(2,2),subplot=1)
    heatmap!(p,x.v./1e3,y.v./1e3,log10.(D.ηs.v'),color=reverse(cgrad(:roma)),
                xlabel="x[km]",ylabel="y[km]",
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                layout=(2,2),subplot=2)
    display(p)
    
    return
    # ------------------------------------------------------------------- #
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    )

    rhs     =   zeros(maximum(Num.Pt))
    χ       =   zeros(maximum(Num.Pt))

    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )

    K       =   Assemblyc(NC, NV, Δ, η₀, VBC, Num)

    rhs     =   updaterhsc( NC, NV, Δ, η₀, D.ρ, g, VBC, Num, rhs )
    
    # Solution --- 
    χ       =   K \ rhs

    # @show size(χ)

    D.vx[:,2:end-1]     .=  χ[Num.Vx]
    D.vy[2:end-1,:]     .=  χ[Num.Vy]
    D.Pt                .=  χ[Num.Pt]

    # Get the velocity on the centroids ---
    for i = 1:NC.x
        for j = 1:NC.y
            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
        end
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

    @show(minimum(D.vc))
    @show(maximum(D.vc))

    quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            quiver=(D.vx[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
            color="white",layout=(3,1),subplot=1)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vxc',
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(3,1),subplot=2)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vyc',
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(3,1),subplot=3)
    display(p)

end

main()