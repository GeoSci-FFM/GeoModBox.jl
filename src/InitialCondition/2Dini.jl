@doc raw""" 
    IniTemperature!(type,M,NC,Δ,D,x,y;Tb=1000,Ta=1200,Ampl=200,σ=0.05)

Function to setup an initial temperature condition for a two dimensional 
problem. The temperature is defined on the centroids of a regular finite 
difference grid. 

    type    : 
    M       : 
    NC      : 
    Δ       : 
    D       : 
    x       : 
    y       : 

    Tb      : Hintergrund Temperatur
    Ta      : Amplitude der Anomalie
    σ       : Breite der Gaussian Anomalie

Possible initial temperature conditions are: 

    1) Circle/elliptical anomaly
    2) Gaussian anomaly
    3) Block anomaly
"""

function IniTemperature!(type,M,NC,Δ,D,x,y;Tb=600,Ta=1200,σ=0.1)

    if type==:circle 
        # Bereich der Anomalie ---       
        ri          =   .2
        xc          =   (M.xmin+M.xmax)/4
        yc          =   (M.ymin+M.ymax)/2
        α           =   0.0
        a_ell       =   .2*(M.ymin+M.ymax)
        b_ell       =   .2*(M.ymin+M.ymax)
        for i = 1:NC.x+2, j = 1:NC.y+2
            x_ell   =  x.cew[i]*cosd(α) + y.cns[j]*sind(α)
            y_ell   =  -x.cew[i]*sind(α) + y.cns[j]*cosd(α)
            Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
            if Elli <= ri 
                D.T_ex[i,j]    =   Ta
            else
                D.T_ex[i,j]    =   Tb
            end
        end    
        D.Tmax[1]   =   maximum(D.T_ex)
        D.Tmin[1]   =   minimum(D.T_ex)
        D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
    elseif type==:gaussian        
        # κ           =   1e-6
        # AnalyticalSolution2D!(D.T, x.c, y.c, 0.0, (T0=Ampl,K=κ,σ=σ))
        for i = 1:NC.x+2, j = 1:NC.y+2
            D.T_ex[i,j]    =   Tb + Ta*exp(-((x.cew[i]/((M.ymin+M.ymax)) - 0.20)^2 + 
                                    (y.cns[j]/((M.ymin+M.ymax)) - 0.5)^2)/σ^2)
        end        
        D.Tmax[1]   =   maximum(D.T_ex)
        D.Tmin[1]   =   minimum(D.T_ex)
        D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
    elseif type==:block        
        # Bereich der Temperatur Anomalie ---
        xTl     =   (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/4 - (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/10
        xTr     =   (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/4 + (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/10
        yTu     =   (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/2 - (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/10
        yTo     =   (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/2 + (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/10
        Ta      =   1200
        D.Tmean[1]  =   (Tb + Ta)/2
        # Anfangstemperatur Verteilung ---
        for i = 1:NC.x+2, j = 1:NC.y+2
            if y.cns[j]>=yTu && y.cns[j] <= yTo && x.cew[i]>=xTl && x.cew[i]<=xTr
                D.T_ex[i,j]    =   Ta
            else
                D.T_ex[i,j]    =   Tb
            end
        end        
        D.Tmax[1]   =   maximum(D.T_ex)
    end
    # Assign temperature to regular field ---
    D.T         .=  D.T_ex[2:end-1,2:end-1]
    return D
end    

@doc raw"""
    IniVelocity!(type,D,NC,Δ,M,x,y)
#
#
# ---
"""

function IniVelocity!(type,D,NV,Δ,M,x,y)
    if type==:RigidBody
        # Rigid Body Rotation ---
        for i = 1:NV.x, j = 1:NV.y+1
            D.vx[i,j]  =    ((y.cns[j]-(M.ymax-M.ymin)/2))/(M.ymax-M.ymin)
        end
        for i = 1:NV.x+1, j = 1:NV.y
            D.vy[i,j]  =   -((x.cew[i]-(M.xmax-M.xmin)/2))/(M.ymax-M.ymin)
        end
        
        Radx        =   zeros(size(D.vx))
        Rady        =   zeros(size(D.vy))

        @. Radx     =   sqrt((x.vx2d-(M.xmax-M.xmin)/2)^2 + (y.vx2d-(M.ymax-M.ymin)/2)^2)
        @. Rady     =   sqrt((x.vy2d-(M.xmax-M.xmin)/2)^2 + (y.vy2d-(M.ymax-M.ymin)/2)^2)

        @. D.vx[Radx>(M.xmax-M.xmin)/2-1*Δ.x]     =   0
        @. D.vy[Rady>(M.xmax-M.xmin)/2-1*Δ.x]     =   0

        @. D.vx     =   D.vx/(100*(60*60*24*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100*(60*60*24*365.25))        # [m/s]
        
    elseif type==:ShearCell
        # Convection Cell with a Shear Deformation --- (REF?!)
        for i = 1:NV.x, j = 1:NV.y+1
            D.vx[i,j]   =   -sin(π*x.v[i])*cos(π*y.cns[j])
        end
        for i = 1:NV.x+1, j = 1:NV.y
            D.vy[i,j]   =   cos(π.*x.cew[i]).*sin(π.*y.v[j])
        end
        @. D.vx     =   D.vx/(100*(60*60*24*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100*(60*60*24*365.25))        # [m/s]
    end
    return D
end