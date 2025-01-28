function ForwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC)
    # Function to solve 2D heat diffusion equation using the explicit finite
    # difference scheme
    # Q - Waermeproduktionsrate pro Volumen [W/m^3]
    # ------------------------------------------------------------------- #
    
    sx      = κ * Δt / Δx^2
    sz      = κ * Δt / Δy^2

    # D.T_ex[2:end-1,2:end-1]     .=  D.T

    # Temperature at the ghost nodes ------------------------------------ #
    # West boundary ---
    D.T_ex[1,2:end-1]   .= (BC.type.W==:Dirichlet) .* (2 .* BC.val.W .- D.T_ex[2,2:end-1]) + 
                            (BC.type.W==:Neumann) .* (D.T_ex[2,2:end-1] .- BC.val.W .* Δx)
    # East boundary ---
    D.T_ex[end,2:end-1] .= (BC.type.E==:Dirichlet) .* (2 .* BC.val.E .- D.T_ex[end-1,2:end-1]) + 
                            (BC.type.E==:Neumann) .* (D.T_ex[end-1,2:end-1] .+ BC.val.E .* Δx)
    # South boundary --- 
    D.T_ex[2:end-1,1]   .= (BC.type.S==:Dirichlet) .* (2 .* BC.val.S .- D.T_ex[2:end-1,2]) + 
                            (BC.type.S==:Neumann) .* (D.T_ex[2:end-1,2] .- BC.val.S .* Δy)
    # Northern boundary ---
    D.T_ex[2:end-1,end] .= (BC.type.N==:Dirichlet) .* (2 .* BC.val.N .- D.T_ex[2:end-1,end-1]) + 
                            (BC.type.N==:Neumann) .* (D.T_ex[2:end-1,end-1] .+ BC.val.N .* Δy)
    # ------------------------------------------------------------------- #
    # Loop over internal nodes ------------------------------------------ #
    for i = 1:NC.x, j = 1:NC.y
        i1 = i+1
        j1 = j+1
        D.T[i,j] = D.T_ex[i1,j1] + 
            sx * (D.T_ex[i1-1,j1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1+1,j1]) + 
            sz * (D.T_ex[i1,j1-1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1,j1+1]) + 
            D.Q[i,j] * Δt / ρ / cp
    end
    # ------------------------------------------------------------------- #
    # Update extended temperature --------------------------------------- #
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #
    
end