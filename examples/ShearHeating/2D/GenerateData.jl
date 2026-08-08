
include("ShearHeatingShearBands.jl")

Diff    =   [
    :dc,
]

θ   =   [
    0,
    0.5,
    1.0
]

Adv     =   [
    :upwind,
    :semilag,
    :tracers,
]

style   =   [
    # :moving,
    :fixed,
    # :max,
]

avg     =   [
    :arithmetic,
    # :geometric,
    # :harmonic,
]

for l in eachindex(avg)
    for k in eachindex(style)
        for i in eachindex(θ)
            for j in eachindex(Adv)
                @show θ[i], Adv[j], style[k], avg[l]
                ShearHeatingShearBands(Diff[1],θ[i],Adv[j],style[k],avg[l])
            end
        end
    end
end

