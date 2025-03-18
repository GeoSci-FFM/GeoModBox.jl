using ExtendableSparse

@doc raw"""
    Assemblyc()
"""
function Assemblyc(NC, NV, Δ, η, BC, Num)

    # Linear system of equation ---
    ndof    =   maximum(Num.Pt)
    K       =   ExtendableSparseMatrix(ndof,ndof)
    dx,dy   =   Δ.x, Δ.y

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NV.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NV.x
            # Pressure ---
            iPE =   Num.Pt[i,j]
            iPW =   Num.Pt[i-1,j]
            # ---
            inS     =   j==1    ? false  : true  
            FSS     =   (j==1    && BC.type.S==:freeslip) ? 1. : 0.
            NSS     =   (j==1    && BC.type.S==:noslip) ? 1. : 0.
            inN     =   j==NC.y ? false  : true   
            FSN     =   (j==NC.y && BC.type.N==:freeslip) ? 1. : 0.
            NSN     =   (j==NC.y && BC.type.N==:noslip) ? 1. : 0.
            if inS K[ii,iS] =   η / dy^2 end
            K[ii,iW]    =   η / dx^2
            K[ii,iC]    =   - 2.0 * η / dx^2 - (2.0 - FSS - FSN + NSS + NSN) * η / dy^2
            K[ii,iE]    =   η / dx^2
            if inN K[ii,iN] =   η / dy^2 end
            K[ii,iPW]   = 1 / dx
            K[ii,iPE]   = -1 / dx
        end
    end
    # ------------------------------------------------------------------- #
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        if j == 1 || j == NV.y
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NC.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NC.x
            # Pressure ---
            iPS =   Num.Pt[i,j-1]   
            iPN =   Num.Pt[i,j]
            # ---
            inW     =   i==1    ? false  : true  
            FSW     =   (i==1    && BC.type.W==:freeslip) ? 1. : 0.
            NSW     =   (i==1    && BC.type.W==:noslip) ? 1. : 0.
            inE     =   i==NC.x ? false  : true   
            FSE     =   (i==NC.x && BC.type.E==:freeslip) ? 1. : 0.
            NSE     =   (i==NC.x && BC.type.E==:noslip) ? 1. : 0.            
            K[ii,iS]    =   η / dy^2
            if inW K[ii,iW] =  η / dx^2 end
            K[ii,iC]    =   - (2.0 - FSW - FSE + NSW + NSE) * η / dx^2 - 2.0 * η / dy^2
            if inE K[ii,iE] = η / dx^2 end
            K[ii,iN]    =   η / dy^2
            K[ii,iPS]   =   1 / dy
            K[ii,iPN]   =   - 1 / dy
        end
    end
    # ------------------------------------------------------------------- #
    # Continuity equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NC.y
        # Equation number ---
        ii  =   Num.Pt[i,j]
        # Stencil ---
        iW  =   Num.Vx[i,j]
        iE  =   Num.Vx[i+1,j]
        iS  =   Num.Vy[i,j]
        iN  =   Num.Vy[i,j+1]
        # Linear system coefficients
        K[ii,iW]    =   -1 / dx
        K[ii,iE]    =   1 / dx
        K[ii,iS]    =   -1 / dy 
        K[ii,iN]    =   1 / dy
    end
    return flush!(K)
end

@doc raw"""
    updaterhsc()
"""
function updaterhsc(NC, NV, Δ, η, ρ, g, BC, Num, rhs)

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        rhs[ii]     =   0.0
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            rhs[ii]     =   0.0
        else            
            # ---
            NSS     =   (j==1    && BC.type.S==:noslip) ? 1. : 0.
            NSN     =   (j==NC.y && BC.type.N==:noslip) ? 1. : 0.            
            # ---
            rhs[ii] +=  -2.0 * η * BC.val.S[i] / Δ.y^2 * NSS -
                        2.0 * η * BC.val.N[i] / Δ.y^2 * NSN
        end
    end
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        rhs[ii]     =   0.0
        if j == 1 || j == NV.y
            # North and South boundary ---
            # Free Slip && No Slip: vy = 0 
            rhs[ii]     =   0.0
        else            
            # ---
            NSW     =   (i==1    && BC.type.W==:noslip) ? 1. : 0.
            NSE     =   (i==NC.x && BC.type.E==:noslip) ? 1. : 0.            
            # ---
            rhs[ii] +=  g * ((ρ[i,j] + ρ[i,j-1]) / 2.0) - 
                        2.0 * η * BC.val.W[i] / Δ.x^2 * NSW -
                        2.0 * η * BC.val.E[i] / Δ.x^2 * NSE
        end
    end
    return rhs
end