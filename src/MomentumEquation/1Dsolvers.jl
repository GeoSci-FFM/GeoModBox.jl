using ExtendableSparse

@doc raw"""
    ComputeStokesResiduals1D! ()
"""
function ComputeStokesResiduals1D!( D, ∂P∂x, Δy, BC)
    @. D.vₓₑ[2:end-1]     =   D.vₓ
    # Bottom ---
    D.vₓₑ[1]            =   (BC.type.S==:Dirichlet)*(2*BC.val.S - D.vₓₑ[2]) + 
                            (BC.type.S==:Neumann)*(D.vₓₑ[2] - Δy*BC.val.S)  
    # Top --- 
    D.vₓₑ[end]          =   (BC.type.N==:Dirichlet)*(2*BC.val.N - D.vₓₑ[end-1]) + 
                            (BC.type.N==:Neumann)*(D.vₓₑ[end-1] - Δy*BC.val.N)  
    # ---
    @. D.∂vₓ∂y          =   (D.vₓₑ[2:end] - D.vₓₑ[1:end-1]) / Δy
    @. D.τxy            =   D.η * D.∂vₓ∂y
    @. D.∂τxy∂y         =   (D.τxy[2:end] - D.τxy[1:end-1]) / Δy
    @. D.R              =   -∂P∂x + D.∂τxy∂y
end

@doc raw"""
    AssembleStokesMatrix1D()
"""
function AssembleStokesMatrix1D(nc, η, Δy, BC, K)
    # Zusammenstellen der Koeffizientenmatrix ------------------------------- #
    for i = 1:nc
        a   =   η[i] / Δy^2.0
        b   =   -(η[i]+η[i+1]) / Δy^2.0
        c   =   η[i+1] / Δy^2.0
        # Gleichungsnummer ---
        ii  =   i
        # Stempel ---
        iN  =   ii + 1      #   Norden
        iC  =   ii          #   Zentral
        iS  =   ii - 1      #   Süden
        # Ränder ---
        # Falls ein Süd Index gebrauch wird ---
        inS    =  i==1    ? false  : true
        DirS   = (i==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (i==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inN    =  i==nc ? false  : true
        DirN   = (i==nc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (i==nc && BC.type.N==:Neumann  ) ? 1. : 0.
        if inS K[ii,iS]    = a end
            K[ii,iC]       = b + (NeuS - DirS)*a + (NeuN - DirN)*c
        if inN K[ii,iN]    = c end    
    end
    return flush!(K)
end

@doc raw"""
    Stokes_1D()
"""
function Stokes_1D_direct(vₓ,η,Δy,nc,BC,K,rhs)

    # Zusammenstellen der Koeffizientenmatrix ------------------------------- #
    for i = 1:nc
        a   =   η[i] / Δy^2.0
        b   =   -(η[i]+η[i+1]) / Δy^2.0
        c   =   η[i+1] / Δy^2.0
        # Gleichungsnummer ---
        ii  =   i
        # Stempel ---
        iN  =   ii + 1      #   Norden
        iC  =   ii          #   Zentral
        iS  =   ii - 1      #   Süden
        # Ränder ---
        # Falls ein Süd Index gebrauch wird ---
        inS    =  i==1    ? false  : true
        DirS   = (i==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (i==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inN    =  i==nc ? false  : true
        DirN   = (i==nc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (i==nc && BC.type.N==:Neumann  ) ? 1. : 0.
        if inS K[ii,iS]    = a end
            K[ii,iC]       = b + (NeuS - DirS)*a + (NeuN - DirN)*c
        if inN K[ii,iN]    = c end    
        # Änderung der rechten Seite ---
        rhs[i]      +=  a*BC.val.S*Δy * NeuS - 
                            2*a*BC.val.S * DirS - 
                            c*BC.val.N*Δy * NeuN - 
                            2*c*BC.val.N * DirN
    end
    # ----------------------------------------------------------------------- #
    # Lösung des Gleichungssystems ------------------------------------------ #
    vₓ      .=   K \ rhs
    # ----------------------------------------------------------------------- #
    return vₓ

end