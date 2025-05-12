using Base.Threads 

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
@views function IniTemperature!(type,M,NC,Δ,D,x,y;Tb=600.0,Ta=1200.0,σ=0.1)
    if type==:circle 
        # Circle shaped anomaly ---
        # Bereich der Anomalie ---       
        ri          =   .2
        xc          =   (M.xmin[1]+M.xmax[1])/4
        yc          =   (M.ymin[1]+M.ymax[1])/2
        α           =   0.0
        a_ell       =   .2*(M.ymin[1]+M.ymax[1])
        b_ell       =   .2*(M.ymin[1]+M.ymax[1])
        @threads for i = 1:NC.x+2 
            for j = 1:NC.y+2
                x_ell   =  x.ce[i]*cosd(α) + y.ce[j]*sind(α)
                y_ell   =  -x.ce[i]*sind(α) + y.ce[j]*cosd(α)
                Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
                if Elli <= ri 
                    D.T_ex[i,j]    =   Ta
                else
                    D.T_ex[i,j]    =   Tb
                end
            end
        end    
    elseif type==:gaussian        
        # κ           =   1e-6
        # AnalyticalSolution2D!(D.T, x.c, y.c, 0.0, (T0=Ampl,K=κ,σ=σ))
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                D.T_ex[i,j]    =   Tb + Ta*exp(-((x.ce[i]/((M.ymin[1]+M.ymax[1])) - 0.20)^2 + 
                                    (y.ce[j]/((M.ymin[1]+M.ymax[1])) - 0.5)^2)/σ^2)
            end
        end        
    elseif type==:block        
        # Bereich der Temperatur Anomalie ---
        xTl     =   M.xmin[1] + (M.xmax[1]-M.xmin[1])/8.0 # (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/4 - (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/10
        xTr     =   xTl + (M.xmax[1]-M.xmin[1])/10.0 # (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/4 + (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/10
        yTu     =   M.ymin[1] + (M.ymax[1]-M.ymin[1])/2.0 - (M.ymax[1]-M.ymin[1])/10.0  # (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/2 - (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/10
        yTo     =   M.ymin[1] + (M.ymax[1]-M.ymin[1])/2.0 + (M.ymax[1]-M.ymin[1])/10.0 #(abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/2 + (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/10
        # xTl     =   (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/4 - (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/10
        # xTr     =   (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/4 + (abs(M.xmin[1]-Δ.x/2)+abs(M.xmax[1]+Δ.x/2))/10
        # yTu     =   (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/2 - (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/10
        # yTo     =   (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/2 + (abs(M.ymin[1]-Δ.y/2)+abs(M.ymax[1]+Δ.y/2))/10
        # Anfangstemperatur Verteilung ---
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                if y.ce[j]>=yTu && y.ce[j] <= yTo && x.ce[i]>=xTl && x.ce[i]<=xTr
                    D.T_ex[i,j]    =   Ta
                else
                    D.T_ex[i,j]    =   Tb
                end
            end
        end        
    elseif type==:linear
        Ttop    =   Ta
        Tbot    =   Tb
        Tgrad   =   (Tbot-Ttop)/(M.ymax[1]-M.ymin[1])         # [ K/m ]
        @show Tgrad
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                D.T_ex[i,j] = -Tgrad*(y.ce[j]) + Ttop
            end
        end
    elseif type==:lineara
        # Bereich der Anomalie ---       
        ri          =   .3
        xc          =   (M.xmin[1]+M.xmax[1])/4
        yc          =   (M.ymin[1]+M.ymax[1])/2
        α           =   0.0
        a_ell       =   .6*(M.ymin[1]+M.ymax[1])
        b_ell       =   .2*(M.ymin[1]+M.ymax[1])
        # Linear with a gaussian anomaly
        Ttop    =   Ta
        Tbot    =   Tb
        Tgrad   =   (Tbot-Ttop)/(M.ymax[1]-M.ymin[1])         # [ K/m ]
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                x_ell   =  x.ce[i]*cosd(α) + y.ce[j]*sind(α)
                y_ell   =  -x.ce[i]*sind(α) + y.ce[j]*cosd(α)
                Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
                D.T_ex[i,j] = -Tgrad*(y.ce[j]) + Ttop
                if Elli <= ri
                    D.T_ex[i,j]    +=   0.2*D.T_ex[i,j]
                end
            end
        end
    end
    D.Tmax[1]   =   maximum(D.T_ex)
    D.Tmin[1]   =   minimum(D.T_ex)
    D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
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
@views function IniVelocity!(type,D,NV,Δ,M,x,y)
    if type==:RigidBody
        # Rigid Body Rotation ---
        # We assume a maximum and minimum velocity of 0.5 cm/a, respectively! 
        @threads for i = 1:NV.x
            for j = 1:NV.y+1
                D.vx[i,j]  =    ((y.ce[j]-(M.ymax[1]-M.ymin[1])/2))/(M.ymax[1]-M.ymin[1])
            end
        end
        @threads for i = 1:NV.x+1
            for j = 1:NV.y
                D.vy[i,j]  =   -((x.ce[i]-(M.xmax[1]-M.xmin[1])/2))/(M.ymax[1]-M.ymin[1])
            end
        end
        
        Radx        =   zeros(size(D.vx))
        Rady        =   zeros(size(D.vy))

        @. Radx     =   sqrt((x.vx2d-(M.xmax[1]-M.xmin[1])/2)^2 + (y.vx2d-(M.ymax[1]-M.ymin[1])/2)^2)
        @. Rady     =   sqrt((x.vy2d-(M.xmax[1]-M.xmin[1])/2)^2 + (y.vy2d-(M.ymax[1]-M.ymin[1])/2)^2)

        @. D.vx[Radx>(M.xmax[1]-M.xmin[1])/2-1*Δ.x]     =   0
        @. D.vy[Rady>(M.xmax[1]-M.xmin[1])/2-1*Δ.x]     =   0

        @. D.vx     =   D.vx/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        
    elseif type==:ShearCell
        # Convection Cell with a Shear Deformation --- (REF?!)
        @threads for i = 1:NV.x 
            for j = 1:NV.y+1
                D.vx[i,j]   =   -sin(π*(x.v[i]/(M.xmax[1]-M.xmin[1])))*
                                    cos(π*y.ce[j]/(M.ymax[1]-M.ymin[1]))
            end
        end
        @threads for i = 1:NV.x+1 
            for j = 1:NV.y
                D.vy[i,j]   =   cos(π*x.ce[i]/(M.xmax[1]-M.xmin[1]))*
                                    sin(π*y.v[j]/(M.ymax[1]-M.ymin[1]))
            end
        end
        @. D.vx     =   D.vx/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
    end
    return D
end

@doc raw"""
    IniPhase!
"""
@views function IniPhase!(type,D,M,x,y,NC;phase=0)

    if type==:block
        # Bereich der Anomalie ---
        xL      =   2/5 * (M.xmax[1]-M.xmin[1])
        xR      =   3/5 * (M.xmax[1]-M.xmin[1])
        yO      =   0.1 * (M.ymin[1]-M.ymax[1])
        yU      =   0.3 * (M.ymin[1]-M.ymax[1])        
        
        # Phase ---
        for i = 1:NC.x
            for j = 1:NC.y
                if y.c[j]>=yU && y.c[j] <= yO && x.c[i]>=xL && x.c[i]<=xR
                    D.p[i,j]    =   phase[2]    #   anomaly 
                else
                    D.p[i,j]    =   phase[1]    #   background
                end
            end
        end
        D.p_ex[2:end-1,2:end-1]     .=  D.p
        D.p_ex[1,:]     .=   D.p_ex[2,:]
        D.p_ex[end,:]   .=   D.p_ex[end-1,:]
        D.p_ex[:,1]     .=   D.p_ex[:,2]
        D.p_ex[:,end]   .=   D.p_ex[:,end-1]
        # # Viscosity ---
        # if size(η,1)==2
        #     for i = 1:NC.x
        #         for j = 1:NC.y
        #             if y.c[j]>=yU && y.c[j] <= yO && x.c[i]>=xL && x.c[i]<=xR
        #                 D.ηs.v[i,j] =   η[2]
        #             else
        #                 D.ηs.v[i,j] =   η[1]
        #             end
        #         end
        #     end 
        #     @. D.ηs.c = 0.25*(D.ηs.v[1:end-1,1:end-1] + 
        #                     D.ηs.v[2:end-0,1:end-1] + 
        #                     D.ηs.v[1:end-1,2:end-0] + 
        #                     D.ηs.v[2:end-0,2:end-0])
        # end
    end
    return D
end