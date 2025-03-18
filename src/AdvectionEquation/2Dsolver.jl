# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
using Interpolations

@doc raw"""
    upwindc2D!()

    Using the velocity on the centroids! 
"""
function upwindc2D!(P,P_ex,vxc,vyc,NC,Δt,Δx,Δy)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. P     =  P_ex[indx,indy] - 
            (vxc>0)*(vxc*Δt/Δx*(P_ex[indx,indy  ] - P_ex[indx-1,indy])) - 
            (vxc<0)*(vxc*Δt/Δx*(P_ex[indx+1,indy] - P_ex[indx,indy  ])) - 
            (vyc>0)*(vyc*Δt/Δy*(P_ex[indx,indy  ] - P_ex[indx,indy-1])) - 
            (vyc<0)*(vyc*Δt/Δy*(P_ex[indx,indy+1] - P_ex[indx,indy  ]))
    # Update extende temperature field ------------------------------- #
    P_ex[indx,indy]     .=  P
    # ---------------------------------------------------------------- #
end

@doc raw"""
    slfc2D!
"""
function slfc2D!(P,P_ex,P_exo,vxc,vyc,NC,Δt,Δx,Δy)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. P  =   P_exo[indx,indy] - 
        vxc*Δt/Δx*(P_ex[indx+1,indy] - P_ex[indx.-1,indy]) - 
        vyc*Δt/Δy*(P_ex[indx,indy+1] - P_ex[indx,indy-1])
    # Update extende temperature field ------------------------------- #
    @. P_exo            =  P_ex
    P_ex[indx,indy]     .=  P
    # ---------------------------------------------------------------- #
end

@doc raw"""
    semilagc2D!()
"""
function semilagc2D!(P,P_ex,vxc,vyc,vxo,vyo,x,y,Δt)
    # mid-point iteration scheme ---
    if isempty(vxo) || isempty(vyo)
        # Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
        # aktuelle Geschwindigkeit auf die alte Geschwindigkeit
        # uebertragen.
        vxo   =   copy(vxc)
        vyo   =   copy(vyc)
    end
    
    # Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --- 
    vxcm        =   copy(vxc)
    vycm        =   copy(vyc)
    @. vxcm     =   0.5*(vxo + vxc)
    @. vycm     =   0.5*(vyo + vyc)
    
    # Initialisierung der Geschwindigkeit fuer die Iteration ---
    vxi     =   copy(vxc)
    vyi     =   copy(vyc)
    xp      =   copy(x.c2d)
    yp      =   copy(y.c2d)
    
    # Iteration ---
    for k = 1:10
        @. xp  = x.c2d - 0.5*Δt*vxi
        @. yp  = y.c2d - 0.5*Δt*vyi

        itp_linearx  =  linear_interpolation((x.c,y.c),vxcm,extrapolation_bc = Line())
        itp_lineary  =  linear_interpolation((x.c,y.c),vycm,extrapolation_bc = Line())

        vxi         .=  itp_linearx.(xp,yp)
        vyi         .=  itp_lineary.(xp,yp)
    end
    @. xp   =   x.c2d - Δt*vxi
    @. yp   =   y.c2d - Δt*vyi        

    itp_cubic   =   cubic_spline_interpolation((x.ce,y.ce),P_ex)
    P           .=  itp_cubic.(xp,yp)
        
    P_ex[2:end-1,2:end-1]   .=  P
end