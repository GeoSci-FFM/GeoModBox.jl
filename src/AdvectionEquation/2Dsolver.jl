# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
using Interpolations

function upwindc2D!(D,NC,T,Δ)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. D.T     =  D.T_ex[indx,indy] - 
            (D.vxc>0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx,indy] - D.T_ex[indx-1,indy])) - 
            (D.vxc<0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx+1,indy] - D.T_ex[indx,indy])) - 
            (D.vyc>0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy] - D.T_ex[indx,indy-1])) - 
            (D.vyc<0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy+1] - D.T_ex[indx,indy]))
    # Update extende temperature field ------------------------------- #
    D.T_ex[indx,indy]  .=  D.T
    # ---------------------------------------------------------------- #
end

function slfc2D!(D,NC,T,Δ)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. D.T  =   D.T_exo[indx,indy] - 
        D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx+1,indy]-D.T_ex[indx.-1,indy]) - 
        D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy+1]-D.T_ex[indx,indy-1])
    # Update extende temperature field ------------------------------- #
    @. D.T_exo         =  D.T_ex
    D.T_ex[indx,indy]  .=  D.T
    # ---------------------------------------------------------------- #
end

function semilagc2D!(D,vxo,vyo,x,y,T)
    # mid-point iteration scheme ---
    if isempty(vxo) || isempty(vyo)
        # Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
        # aktuelle Geschwindigkeit auf die alte Geschwindigkeit
        # uebertragen.
        vxo   =   copy(D.vxc)
        vyo   =   copy(D.vyc)
    end
    
    # Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --- 
    vxcm        =   copy(D.vxc)
    vycm        =   copy(D.vyc)
    @. vxcm     =   0.5*(vxo + D.vxc)
    @. vycm     =   0.5*(vyo + D.vyc)
    
    # Initialisierung der Geschwindigkeit fuer die Iteration ---
    vxi     =   copy(D.vxc)
    vyi     =   copy(D.vyc)
    xp      =   copy(x.c2d)
    yp      =   copy(y.c2d)
    
    # Iteration ---
    for k = 1:10
        @. xp  = x.c2d - 0.5*T.Δ[1]*vxi
        @. yp  = y.c2d - 0.5*T.Δ[1]*vyi

        itp_linearx  =  linear_interpolation((x.c,y.c),vxcm,extrapolation_bc = Line())
        itp_lineary  =  linear_interpolation((x.c,y.c),vycm,extrapolation_bc = Line())

        vxi         .=  itp_linearx.(xp,yp)
        vyi         .=  itp_lineary.(xp,yp)
    end
    @. xp   =   x.c2d - T.Δ[1]*vxi
    @. yp   =   y.c2d - T.Δ[1]*vyi        

    itp_cubic   =   cubic_spline_interpolation((x.cew,y.cns),D.T_ex)
    D.T         .=  itp_cubic.(xp,yp)
        
    D.T_ex[2:end-1,2:end-1]  .=  D.T
    return D
end