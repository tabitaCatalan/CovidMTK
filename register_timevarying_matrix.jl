#=
Minimal example of working with a registered time varying matrix
=#
using ModelingToolkit, OrdinaryDiffEq

interactions_in_time = rand(2,2,10)

function interaction_matrix(t, i, j) 
    interactions_in_time[i,j, floor(Int, t) + 1] 
end


@variables i j
@register interaction_matrix(t, i, j)
#@register interaction_matrix(t)::Array{Float64,2}
@variables t
@variables X[1:2](t) IM[1:2,1:2](t)

D = Differential(t)


eqs = [
    vec([IM[i,j] ~ interaction_matrix(t, i, j) for i in 1:2, j in 1:2]);
    collect(D.(X) .~ (IM * X));
]


@named system = ODESystem(eqs) 

u0 = [ X[1] => 1., X[2] => 2.]
p = []
prob = ODEProblem(structural_simplify(system), u0, (0.0,9.0), p)
sol = solve(prob)