function CNA2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K1, K2, Num)
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
# ----------------------------------------------------------------------- #

# Define coefficients ---
a       =   κ / 2 / Δx^2
b       =   κ / 2 / Δy^2
c       =   1 / Δt

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
    inW    =  i==1    ? false  : true
    DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
    NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
    # If an East index is required ---
    inE    =  i==NC.x ? false  : true
    DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
    NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
    # If an South index is required
    inS    =  j==1      ? false  : true
    DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
    NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
    # If an North index is required 
    inN    =  j==NC.y   ? false  : true
    DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
    NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
    if inS
        K1[ii,iS]   =   - b
        K2[ii,iS]   =   b
    end
    if inE
        K1[ii,iE]   =   - a
        K2[ii,iE]   =   a            
    end
    K1[ii,iC]       =   c + (2 + DirW + DirE - NeuE - NeuW) * a + (2 + DirS + DirN - NeuS - NeuN) * b
    K2[ii,iC]       =   c - (2 + DirW + DirE - NeuE - NeuW) * a - (2 + DirS + DirN - NeuS - NeuN) * b
    if inW 
        K1[ii,iW]   =   - a
        K2[ii,iW]   =   a
    end           
    if inN
        K1[ii,iN]   =   - b
        K2[ii,iN]   =   b
    end
end
# ------------------------------------------------------------------- #
# Berechnung der rechten Seite -------------------------------------- #
rhs     .=   K2 * reshape(D.T,NC.x*NC.y) .+ reshape(D.Q,NC.x*NC.y)./ρ./cp
# ------------------------------------------------------------------- #        
# Aenderung der rechten Seite durch die Randbedingungen ------------- #    
for i = 1:NC.x, j = 1:NC.y     
    ii      =   Num.T[i,j]
    # Boundaries         
    DirW    =   (i==1       && BC.type.W==:Dirichlet) ? 1. : 0.
    NeuW    =   (i==1       && BC.type.W==:Neumann  ) ? 1. : 0.
    DirE    =   (i==NC.x    && BC.type.E==:Dirichlet) ? 1. : 0.
    NeuE    =   (i==NC.x    && BC.type.E==:Neumann  ) ? 1. : 0.
    DirS    =   (j==1       && BC.type.S==:Dirichlet) ? 1. : 0.
    NeuS    =   (j==1       && BC.type.S==:Neumann  ) ? 1. : 0.
    DirN    =   (j==NC.y    && BC.type.N==:Dirichlet) ? 1. : 0.
    NeuN    =   (j==NC.y    && BC.type.N==:Neumann  ) ? 1. : 0.
    
    # Update right hand side ---
    rhs[ii]     = rhs[ii] +
                    4*a*BC.val.W[j] * DirW + 
                    4*a*BC.val.E[j] * DirE +
                    4*b*BC.val.S[i] * DirS +
                    4*b*BC.val.N[i] * DirN - 
                    2*a*BC.val.W[j]*Δx * NeuW +
                    2*a*BC.val.E[j]*Δx * NeuE -
                    2*b*BC.val.S[j]*Δy * NeuS + 
                    2*b*BC.val.N[j]*Δy * NeuN
end
# ------------------------------------------------------------------- #    
# Compute new temperature ------------------------------------------- #
D.T[:]      .=    K1 \ rhs[:]
# ------------------------------------------------------------------- #

end