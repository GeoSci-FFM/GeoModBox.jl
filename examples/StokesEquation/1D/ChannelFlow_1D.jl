using Plots, ExtendableSparse, Printf, LinearAlgebra
using GeoModBox.MomentumEquation.OneD

function ChannelFlow_1D()
# Model Parameter ------------------------------------------------------- #
M   =   (
    ymin        =   -400.0e3,           #   Depth [ m ]
    ymax        =   0.0e3,              
)
I   =   (
    vₓ₀         =   5/100/31536000,     #   Velocity top [ m/s ]
    η₀          =   1.0e21,             #   Viscosity top [ Pa s ]
    η₁          =   1.0e18,             #   Viscosity bottom [ Pa s ]
    ∂P∂x        =   -0.1e1,             #   horizontal pressure gradient [ Pa/m ]
)
I1  =   (
    m           =   I.η₁ / I.η₀,        #   Viscosity ratio
)
I   =   merge(I,I1)
# ----------------------------------------------------------------------- #
# Numerical Parameter --------------------------------------------------- #
NC  =   (
    y   =   100,                        #   Number of centroids
)
NV  =   (
    y   =   NC.y + 1,                   #   Number of vertices
)
Δ   =   (
    y   =   (M.ymax-M.ymin)/NC.y,       #   Grid resolution
)
# Iterations ---
niter   =   10
ϵ       =   1e-10
# ----------------------------------------------------------------------- #
# Allocations ----------------------------------------------------------- #
D   =   (
    η       =   zeros(NC.y+1),
    vₓ      =   zeros(NC...),
    vₓₐ     =   zeros(NC...),
    Δvₓ     =   zeros(NC...),
    vₓₑ     =   zeros(NC.y+2),
    δvₓ     =   zeros(NC...),
    R       =   zeros(NC...),
    ∂τxy∂y  =   zeros(NC...),
    τxy     =   zeros(NV...),
    ∂vₓ∂y   =   zeros(NV...),
)
# ----------------------------------------------------------------------- #
# Grid ------------------------------------------------------------------ #
y   =   (
    c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
    v   =   LinRange(M.ymin,M.ymax,NV.y),
)
@. D.η  =   I.η₀ * exp(-log(I.m)* y.v / (M.ymax-M.ymin))
# ----------------------------------------------------------------------- #
# Analytical Solution --------------------------------------------------- #
if I.m  == 1.0
    @. D.vₓₐ    =   1.0/2.0/I.η₀ * I.∂P∂x * 
                    (y.c^2 + (M.ymax-M.ymin).*y.c) + 
                    I.vₓ₀*y.c/(M.ymax-M.ymin) + I.vₓ₀
else
    @. D.vₓₐ    =   -I.∂P∂x * (M.ymax-M.ymin) / I.η₀ / log(I.m) / (I.m-1) * 
        (-y.c * (I.m^((y.c + (M.ymax-M.ymin))/(M.ymax-M.ymin)) - 
        I.m^(y.c/(M.ymax-M.ymin))) + (M.ymax-M.ymin)*(I.m^(y.c/(M.ymax-M.ymin)) - 1)) + 
        I.vₓ₀ / (I.m-1) * (I.m ^ ((y.c+(M.ymax-M.ymin))/(M.ymax-M.ymin)) - 1)
end
# ----------------------------------------------------------------------- #
# Boundary Conditions --------------------------------------------------- #
VBC     =   (
    type    = (S=:Dirichlet, N=:Dirichlet),
    val     = (S=0.0,N=I.vₓ₀)
)
# ----------------------------------------------------------------------- #
# Coeffficientmatrix ---------------------------------------------------- #
Num     =   (Vx=1:NC.y,)
ndof    =   maximum(Num.Vx)
K       =   ExtendableSparseMatrix(ndof,ndof)
# ----------------------------------------------------------------------- #
# Solution -------------------------------------------------------------- #
for iter=1:niter
    # Evaluate residual ---
    ComputeStokesResiduals1D!(D, I.∂P∂x, Δ.y, VBC)
    @printf("||R|| = %1.4e\n", norm(D.R)/length(D.R))
    norm(D.R)/length(D.R) < ϵ ? break : nothing
    # Assemble linear system ---
    K  = AssembleStokesMatrix1D(NC.y, D.η, Δ.y, VBC, K)
    # Solve for temperature correction: Back substitutions ---
    D.δvₓ .= -(K\D.R[:]) 
    # Update temperature ---
    @. D.vₓ += D.δvₓ[Num.Vx]
end
# ----------------------------------------------------------------------- #
# Deviation from the analytical solution -------------------------------- #
@. D.Δvₓ    =   ((D.vₓₐ - D.vₓ) / D.vₓₐ) * 100.0
# ----------------------------------------------------------------------- #
# Plotting -------------------------------------------------------------- #
q   =   plot(D.vₓ,y.c./1e3,label="numerical",
            xlabel="vₓ", ylabel="y [km]",
            title="Velocity Profile",
            xlim=(0,I.vₓ₀*1.5),ylim=(M.ymin/1e3,M.ymax/1e3),
            layout=(1,3),subplot=1)
plot!(q,D.vₓₐ,y.c./1e3,label="analytical",linestyle=:dash,
        subplot=1)
plot!(q,D.Δvₓ,y.c./1e3,label="",
        xlabel="Δvₓ [%]",ylabel="z [km]",
        title="Error",
        subplot=2)
plot!(q,log10.(D.η),y.v./1e3,label="",
        xlabel="η [Pa s]",ylabel="z [km]",
        title="Viscosity",
        subplot=3)
display(q)
# ----------------------------------------------------------------------- #
savefig("./examples/StokesEquation/1D/Results/ChannelFlow.png")
end

ChannelFlow_1D()