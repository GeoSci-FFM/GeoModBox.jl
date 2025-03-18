using Base.Threads
using Statistics, Printf

mutable struct TMarkers
    x       ::  Array{Float64,1}
    y       ::  Array{Float64,1}
    T       ::  Array{Float64,1}
    phase   ::  Array{Int64,1}
end

mutable struct Markers
    x       ::  Array{Float64,1}
    y       ::  Array{Float64,1}
    phase   ::  Array{Int64,1}
end

@doc raw"""
    IniTracer2D(nmx,nmy,Δ,M,NC,noise)
"""
@views function IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,ini,phase)
    
    nmark   =   nmx*nmy*NC.x*NC.y

    # Initialise markers ---
    Δxm, Δym  =   Δ.x/nmx, Δ.y/nmy 
    
    xm  =   LinRange(M.xmin+Δxm/2.0, M.xmax-Δxm/2.0, NC.x*nmx)
    ym  =   LinRange(M.ymin+Δym/2.0, M.ymax-Δym/2.0, NC.y*nmy)        
    
    (xmi,ymi) = ([x for x=xm,y=ym], [y for x=xm,y=ym])
    
    # Over allocate markers ---
    if Aparam==:thermal        
        Ma  =   TMarkers( vec(xmi), vec(ymi), zeros(Float64, nmark), zeros(Int64, nmark) )
    elseif Aparam==:phase
        Ma  =   Markers( vec(xmi), vec(ymi), zeros(Int64, nmark) )
    end

    # add noise ---
    if noise==1
        @threads for k=1:nmark
            Ma.x[k] += (rand()-0.5)*Δxm
            Ma.y[k] += (rand()-0.5)*Δym
        end
    end

    if ini==:block
        # Geometry of the block anomaly ---
        xL      =   2/5 * (M.xmax-M.xmin)
        xR      =   3/5 * (M.xmax-M.xmin)
        yO      =   0.1 * (M.ymin-M.ymax)
        yU      =   0.3 * (M.ymin-M.ymax)        
        
        # phase ---
        for k = 1:nmark
            if ymi[k]>=yU && ymi[k]<=yO && xmi[k]>=xL && xmi[k]<=xR
                Ma.phase[k]     =   phase[2]    #   anomaly 
            else
                Ma.phase[k]     =   phase[1]    #   background
            end
        end
    end
    return Ma
end

@doc raw"""
    CountMPC()
"""
@views function CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,it)
     # Disable markers outside of the domain
     @threads for k=1:nmark
        if (Ma.x[k]<M.xmin || Ma.x[k]>M.xmax || Ma.y[k]<M.ymin || Ma.y[k]>M.ymax) 
            @inbounds Ma.phase[k] = -1
        end
    end
    # How many are outside? save indices for reuse  
    nmark_out_th    =   zeros(Int64, nthreads())  
    @threads for k=1:nmark
        if Ma.phase[k] == -1
            nmark_out_th[threadid()] += 1
        end
    end
    nmark_out       =   0
    for ith=1:nthreads()
        nmark_out += nmark_out_th[ith]
    end
    @printf("%d markers out\n", nmark_out[1])

    # Initialize marker per cell per thread array ---
    @threads for j = 1:NC.y
        for i = 1:NC.x
            for ith=1:nthreads()
                MPC.th[ith,i,j] = 0.0
            end
        end
    end

    # Count marker per cell per thread ---
    @threads for k=1:nmark
        if (Ma.phase[k]>=0)
            # Get the column:
            dstx = Ma.x[k] - x.c[1]
            i    = Int64(round(ceil( (dstx/Δ.x) + 0.5)))
            # Get the line:
            dsty = Ma.y[k] - y.c[1]
            j    = Int64(round(ceil( (dsty/Δ.y) + 0.5)))
            # Increment cell count
            MPC.th[threadid(),i,j] += 1.0
        end
    end

    @threads for j=1:NC.y
        for i=1:NC.x
            for ith=1:nthreads()
                if ith == 1 
                    MPC.c[i,j] = 0.0
                end
                MPC.c[i,j] += MPC.th[ith,i,j]
            end
        end
    end

    MPC.min[it]         = minimum(MPC.c)
    MPC.max[it]         = maximum(MPC.c)
    MPC.mean[it]        = mean(MPC.c)
    #MPC.tot_reseed[it] = nmark_add
    return Ma
end

@doc raw"""
    VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)
"""
@views function VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)
    # Interpolate vx
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] -  x.v[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x
        i = NC.x
    end
    if j<1
        j = 1
    elseif j> NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] -  x.v[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute vx velocity for the top and bottom of the cell --------- #
    vxm13   =   Vx[i,j  ] * (1-Δxmj/Δ.x) + Vx[i+1,j  ]*Δxmj/Δ.x
    vxm24   =   Vx[i,j+1] * (1-Δxmj/Δ.x) + Vx[i+1,j+1]*Δxmj/Δ.x
    # ---------------------------------------------------------------- #
    if new==1 
        if dxmj/dx>=0.5
            if i<ncx
                vxm13 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i,j  ] - 2.0*Vx[i+1,j  ] + Vx[i+2,j  ])
                vxm24 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i,j+1] - 2.0*Vx[i+1,j+1] + Vx[i+2,j+1])
            end
        else
            if i>1
                vxm13 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i-1,j  ] - 2.0*Vx[i,j  ] + Vx[i+1,j  ])
                vxm24 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i-1,j+1] - 2.0*Vx[i,j+1] + Vx[i+1,j+1])
            end
        end
    end
    # ----------------------------------------------------------------- #
    # Compute vx ------------------------------------------------------ #
    vxm     =   (1-Δymi/Δ.y) * vxm13 + (Δymi/Δ.y) * vxm24
    # ----------------------------------------------------------------- #
    return vxm
end

@doc raw"""
    VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)
"""
@views function VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)
    # Interpolate vy
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] -  y.v[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j = 1
    elseif j>NC.y
        j = NC.y
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    dxmj = Ma.x[k] - x.ce[i]
    dymi = Ma.y[k] -  y.v[j]
    # ---------------------------------------------------------------- #
    # Compute vy velocity for the left and right of the cell --------- #
    vym12 = Vy[i,j  ]*(1-dymi/Δ.y) + Vy[i  ,j+1]*dymi/Δ.y
    vym34 = Vy[i+1,j]*(1-dymi/Δ.y) + Vy[i+1,j+1]*dymi/Δ.y
    # ---------------------------------------------------------------- #
    if new==1 
        if dymi/dy>=0.5
            if j<ncy
                vym12 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i,j  ] - 2.0*Vy[i,j+1  ] + Vy[i,j+2  ]);
                vym34 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i+1,j] - 2.0*Vy[i+1,j+1] + Vy[i+1,j+2]);
            end      
        else
            if j>1
                vym12 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i,j-1  ] - 2.0*Vy[i,j  ] + Vy[i,j+1  ]);
                vym34 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i+1,j-1] - 2.0*Vy[i+1,j] + Vy[i+1,j+1]);
            end
        end
    end
    # ---------------------------------------------------------------- #
    # Compute vy ----------------------------------------------------- #
    vym = (1-dxmj/Δ.x)*vym12 + (dxmj/Δ.x)*vym34
    # ---------------------------------------------------------------- #
    return vym
end

@doc raw"""
    VxVyFromPrNodes(Vxp ,Vyp, k, Ma, x, y, Δ, NC )
"""
@views function VxVyFromPrNodes(Vxp ,Vyp, k, Ma, x, y, Δ, NC)
    # Interpolate vx, vy
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j=1
    elseif j>NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute vx, vy velocity ---------------------------------------- #
    vxm = Vxp[i,j]*wtmij + Vxp[i,j+1]*wtmi1j + Vxp[i+1,j]*wtmij1 + Vxp[i+1,j+1]*wtmi1j1
    vym = Vyp[i,j]*wtmij + Vyp[i,j+1]*wtmi1j + Vyp[i+1,j]*wtmij1 + Vyp[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return vxm, vym
end

@doc raw"""
    FromCtoM(Prop, Ma, x, y, Δ, NC)
"""
@views function FromCtoM(Prop, k, Ma, x, y, Δ, NC)
    # Interpolate Property from Centroids to Marker ---
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j=1
    elseif j>NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute Marker Property ---------------------------------------- #
    Propm   =   Prop[i,j]*wtmij + 
                    Prop[i,j+1]*wtmi1j + 
                    Prop[i+1,j]*wtmij1 + 
                    Prop[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return Propm
end

@doc raw"""
    Markers2Cells(Ma,nmark,PG,weight,x,y,Δ,param)
"""
@views function Markers2Cells(Ma,nmark,PG_th,PG,weight_th,weight,x,y,Δ,param,param2)
    PG      .*=     0.0
    weight  .*=     0.0
    if param==:thermal
        PM  =       copy(Ma.T)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PG_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Get the column:
                    dstx = Ma.x[k] - x.c[1]
                    i = ceil(Int, dstx / Δ.x + 0.5)
                    # Get the line:
                    dsty = Ma.y[k] - y.c[1]
                    j = ceil(Int, dsty /  Δ.y + 0.5)
                    # Relative distances
                    Δxm = 2.0 * abs(x.c[i] - Ma.x[k])
                    Δym = 2.0 * abs(y.c[j] - Ma.y[k])
                    # Increment cell counts
                    area = (1.0 - Δxm / Δ.x) * (1.0 - Δym / Δ.y)
                    PG_th[tid][i, j] += PM[k] * area
                    weight_th[tid][i, j] += area
                end
            end
        end
    elseif param==:phase
        PM  =       copy(Ma.phase)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PG_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Get the column:
                    dstx = Ma.x[k] - x.c[1]
                    i = ceil(Int, dstx / Δ.x + 0.5)
                    # Get the line:
                    dsty = Ma.y[k] - y.c[1]
                    j = ceil(Int, dsty /  Δ.y + 0.5)
                    # Relative distances
                    Δxm = 2.0 * abs(x.c[i] - Ma.x[k])
                    Δym = 2.0 * abs(y.c[j] - Ma.y[k])
                    # Increment cell counts
                    area = (1.0 - Δxm / Δ.x) * (1.0 - Δym / Δ.y)
                    PG_th[tid][i, j] += param2[PM[k]+1] * area
                    weight_th[tid][i, j] += area
                end
            end
        end
    end
    
    PG      .= reduce(+, PG_th)
    #phase  .= reduce(+, phase_th)
    weight  .= reduce(+, weight_th)
    PG ./= weight

    return
end

@doc raw"""
    FromCtoM(Prop, Ma, x, y, Δ, NC)
"""
@views function AdvectTracer2D(Ma,nmark,D,x,y,dt,Δ,NC,rkw,rkv,style)
    @threads for k = 1:nmark
        if (Ma.phase[k]>=0)
            x0  =   Ma.x[k]
            y0  =   Ma.y[k]
            vx  =   0.0
            vy  =   0.0
            # Runge-Kutta loop ---
            for rk=1:4
                # Interp velocity from grid ---
                if style == 1 # Bilinear velocity interp (original is Markers_divergence_ALLSCHEMES_RK4.m)
                    vxm = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 0)
                    vym = VyFromVyNodes(D.vy, k, Ma, x, y, Δ, NC, 0)
                elseif style == 2
                    vxx = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 0)
                    vyy = VyFromVxNodes(D.vy, k, Ma, x, y, Δ, NC, 0)
                    vxp, vyp = VxVyFromPrNodes(D.vxc ,D.vyc, k, Ma, x, y, Δ, NC)
                    vxm = itpw*vxp + (1.0-itpw)*vxx
                    vym = itpw*vyp + (1.0-itpw)*vyy
                elseif style == 3
                    vxm = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 1)
                    vym = VyFromVxNodes(D.vy, k, Ma, x, y, Δ, NC, 1)
                end
                # Temporary RK advection steps ---
                Ma.x[k]     =   x0 + rkv[rk]*dt*vxm
                Ma.y[k]     =   y0 + rkv[rk]*dt*vym
                # Average final velocity ---
                vx    += rkw[rk]*vxm
                vy    += rkw[rk]*vym
            end
            # Advect points ---
            Ma.x[k]     =   x0 + rkv[4]*dt*vx
            Ma.y[k]     =   y0 + rkv[4]*dt*vy
        end
    end
    return Ma
end