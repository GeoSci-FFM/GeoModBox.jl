using ExtendableSparse

function ADI2Dc!(T, κ, Δx, Δy, Δt, ρ, cp, NC, BC)
    # Function to solve 2D heat diffusion equation using the alternating direct
    # implicit finite difference scheme.
    # assuming constant k, rho, cp
    # dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
    # ----------------------------------------------------------------------- #
    # Erstellung der durchlaufenden Indizes ----------------------------- #
    # Gleichungssystem fuer ADI Solver:
    Num     = (Th = reshape(1:NC.x*NC.y,NC.x,NC.y),
                Tv = reshape(1:NC.x*NC.y,NC.y,NC.x)')
    
    # Linear System of Equations -------------------------------------------- #
    ndof    =   maximum(Num.Th)
    A       =   ExtendableSparseMatrix(ndof,ndof)
    B       =   ExtendableSparseMatrix(ndof,ndof)
    C       =   ExtendableSparseMatrix(ndof,ndof)
    D       =   ExtendableSparseMatrix(ndof,ndof)
    rhs     =   zeros(ndof)

    # Setup coefficient matrices ---------------------------------------- #
    a       =   κ / Δx^2
    b       =   κ / Δy^2
    c       =   1.0 / (Δt / 2.0)

    # First ADI step, vertical running scheme ---
    # Erster ADI Schritt: A*T^(l+1/2) = B*T^l -> vertical running scheme
    for i = 1:NC.x,j=1:NC.y
        # Equation number ---
        ii      =   Num.Tv[i,j]
        # println(ii)
        # Stencil ---
        iW      =   ii - NC.y
        iS      =   ii - 1
        iC      =   ii
        iN      =   ii + 1
        iE      =   ii + NC.y
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
        if inW
            B[ii,iW]   =  a 
        end
        if inN
            A[ii,iN]   =   - b
        end
        A[ii,iC]       =   c + (2 + DirS + DirN - NeuS - NeuN) * b
        B[ii,iC]       =   c - (2 + DirW + DirE - NeuW - NeuE) * a
        if inS
            A[ii,iS]   =   - b
        end
        if inE
            B[ii,iE]   =   a
        end        
    end

    rhs[:]  .=   B * T.T[:] .+ T.Q[:]./ρ./cp
    
    # Update rhs from the boundary conditions ---
    for i = 1:NC.x, j = 1:NC.y
        ii  =   Num.Tv[i,j]
        # If an West index is required ---
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
        # If an South index is required
        DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an North index is required 
        DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
    
        rhs[ii]     =   rhs[ii] + 
                        2*a*BC.val.W[j] * DirW + 
                        2*a*BC.val.E[j] * DirE + 
                        2*b*BC.val.S[i] * DirS + 
                        2*b*BC.val.N[i] * DirN -
                        a*BC.val.W[j] * NeuW + 
                        a*BC.val.E[j] * NeuE - 
                        b*BC.val.S[i] * NeuS + 
                        b*BC.val.N[i] * NeuN    
    end

    # Temperature at Δt/2 ---
    T.T[:]  .=  A \ rhs[:]

    # Second ADI step, horizontal running scheme ---
    # Zweiter ADI Schritt: C*T^(l+1) = D*T^(l+1/2) -> horizontal running scheme
    for j = 1:NC.y,i=1:NC.x
        # Equation number ---
        ii      =   Num.Th[i,j]
        # println(ii)
        # Stencil ---
        iS      =   ii - NC.x
        iW      =   ii - 1        
        iC      =   ii
        iE      =   ii + 1
        iN      =   ii + NC.x
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
        
        if inN
            D[ii,iN]   =  b
        end
        if inW
            C[ii,iW]   =  - a
        end
        C[ii,iC]       =   c + (2 + DirW + DirE - NeuW - NeuE) * a
        D[ii,iC]       =   c - (2 + DirS + DirN - NeuS - NeuN) * b
        if inE
            C[ii,iE]   =   - a
        end   
        if inS
            D[ii,iS]   =   b
        end        
    end

    # Update rhs to T^{n+1/2} --- 
    rhs[:]  .=   D * T.T[:] .+ T.Q[:]./ρ./cp
    
    # Update rhs from the boundary conditions ---
    for i = 1:NC.x, j = 1:NC.y
        ii  =   Num.Tv[i,j]
        # If an West index is required ---
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
        # If an South index is required
        DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an North index is required 
        DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
    
        rhs[ii]     =   rhs[ii] + 
                        2*a*BC.val.W[j] * DirW + 
                        2*a*BC.val.E[j] * DirE + 
                        2*b*BC.val.S[i] * DirS + 
                        2*b*BC.val.N[i] * DirN -
                        a*BC.val.W[j] * NeuW + 
                        a*BC.val.E[j] * NeuE - 
                        b*BC.val.S[i] * NeuS + 
                        b*BC.val.N[i] * NeuN    
    end

    # Temperature at Δt ---
    T.T[:]  .=  C \ rhs[:]

end