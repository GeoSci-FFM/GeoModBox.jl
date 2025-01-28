using ExtendableSparse

function ComputeResiduals2D!(R, T, T_ex, T0, ∂T, q, ρ, Cp, k, BC, Δ, Δt)
    @. T_ex[2:end-1,2:end-1] = T 
    @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
    @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
    @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
    @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] - Δ.y/k.y[:,end]*BC.val.N)
    @. ∂T.∂x = (T_ex[2:end,2:end-1] - T_ex[1:end-1,2:end-1])/Δ.x
    @. ∂T.∂y = (T_ex[2:end-1,2:end] - T_ex[2:end-1,1:end-1])/Δ.y
    @. q.x   = -k.x * ∂T.∂x
    @. q.y   = -k.y * ∂T.∂y
    @. R     = ρ*Cp*(T - T0)/Δt + (q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y  
end

function AssembleMatrix2D(rho, cp, k, BC, Num, nc, Δ, Δt)
    # Linear system of equation
    ndof   = maximum(Num.T)
    K      = ExtendableSparseMatrix(ndof, ndof)
    dx, dy = Δ.x, Δ.y
    #############################
    #       Heat equation       #
    #############################
    for i=1:nc.x, j=1:nc.y
        # Equation number
        ii = Num.T[i,j]
        # Stencil
        iS = ii - nc.x
        iW = ii - 1
        iC = ii
        iE = ii + 1
        iN = ii + nc.x
        # Boundaries
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        inE    =  i==nc.x ? false  : true   
        DirE   = (i==nc.x && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc.x && BC.type.E==:Neumann  ) ? 1. : 0.
        inS    =  j==1    ? false  : true  
        DirS   = (j==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        inN    =  j==nc.y ? false  : true   
        DirN   = (j==nc.y && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==nc.y && BC.type.N==:Neumann  ) ? 1. : 0.
        # Material coefficient
        kW = k.x[i,j]
        kE = k.x[i+1,j]
        kS = k.y[i,j]
        kN = k.y[i,j+1]
        ρ  = rho[i,j]
        Cp = cp[i,j]
        # Linear system coefficients
        if inS K[ii,iS] = kS .* (DirS + NeuS - 1) ./ dy .^ 2 end
        if inW K[ii,iW] = kW .* (DirW + NeuW - 1) ./ dx .^ 2 end
        K[ii,iC] = Cp .* ρ ./ Δt + (-kN .* (-DirN + NeuN - 1) ./ dy + kS .* (DirS - NeuS + 1) ./ dy) ./ dy + (-kE .* (-DirE + NeuE - 1) ./ dx + kW .* (DirW - NeuW + 1) ./ dx) ./ dx
        if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ dx .^ 2 end
        if inN K[ii,iN] = -kN .* (-DirN - NeuN + 1) ./ dy .^ 2 end
    end
    return flush!(K)
end

 function BackwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K, Num)
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
# ----------------------------------------------------------------------- #
# Define coefficients ---
a   =   κ / Δx^2
b   =   κ / Δy^2
c   =   1 / Δt
# Multiply rhs with 1/Δt and add Q/ρ/cp ---    
rhs  .= reshape(D.T,NC.x*NC.y).*c .+ reshape(D.Q,NC.x*NC.y)./ρ./cp

# Loop over the grid points ---
for i = 1:NC.x, j = 1:NC.y
    # Equation number ---
    ii          =   Num.T[i,j]
    # Stencil ---
    iS          =   ii - NC.x   # South
    iW          =   ii - 1      # West
    iC          =   ii          # Central
    iE          =   ii + 1      # East
    iN          =   ii + NC.x   # North
    # Boundaries ---
    # If an West index is required ---
    inW    =  i==1      ? false  : true
    DirW   = (i==1      && BC.type.W==:Dirichlet) ? 1. : 0.
    NeuW   = (i==1      && BC.type.W==:Neumann  ) ? 1. : 0.
    # If an East index is required ---
    inE    =  i==NC.x   ? false  : true
    DirE   = (i==NC.x   && BC.type.E==:Dirichlet) ? 1. : 0.
    NeuE   = (i==NC.x   && BC.type.E==:Neumann  ) ? 1. : 0.
    # If an South index is required
    inS    =  j==1      ? false  : true
    DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
    NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
    # If an North index is required 
    inN    =  j==NC.y   ? false  : true
    DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
    NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
    # Stencil ---
    if inS K[ii,iS]     = - b end
    if inW K[ii,iW]     = - a end
    K[ii,iC]            =   (2 + DirW + DirE - NeuW - NeuE)*a + (2 + DirS + DirN - NeuS - NeuN) *b + c
    if inE K[ii,iE]     = - a end    
    if inN K[ii,iN]     = - b end
    # Modify right hand side due to boundary conditions ------------- #
    rhs[ii]     = rhs[ii] + 
                    2*a*BC.val.W[j] * DirW +
                    2*a*BC.val.E[j] * DirE +
                    2*b*BC.val.S[i] * DirS +
                    2*b*BC.val.N[i] * DirN -
                    a*BC.val.W[j]*Δx * NeuW  + 
                    a*BC.val.E[j]*Δx * NeuE  - 
                    b*BC.val.S[i]*Δy * NeuS  + 
                    b*BC.val.N[i]*Δy * NeuN 
end
# ------------------------------------------------------------------- #    
# Calculate temperature at new time step ---------------------------- #
D.T[:]  .=   K \ rhs[:]
# ------------------------------------------------------------------- #

end