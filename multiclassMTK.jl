using ModelingToolkit, OrdinaryDiffEq
using Plots: plot, plot!


n = 2 #length(comunas) # number of classes 
m = 2 # number of environments 

@parameters t γₑ β[1:m] N[1:n]
@variables α[1:n](t) S[1:n](t) E[1:n](t) λ[1:n](t) R[1:n](t)


D = Differential(t) 

# using Random
# variations1 = rand(20)/5;
# variations2 = rand(20)/10; 

# function Pmatrix(t)
#     if t >= 20. 
#         index = 20 
#     else 
#         index = floor(Int, t) + 1
#     end 
#     [
#         (0.8 - variations1[index]) (0.2 + variations1[index]);
#         (0.9 - variations2[index]) (0.1 + variations2[index])
#     ]
# end 
# @register Pmatrix(t)

#@register residence_times_matrix(t)
@register residence_times_matrix2(t)

@variables TRM[1:n, 1:m](t)

eqs = [
    [TRM[i] ~ residence_times_matrix2(t)[i] for i in 1:n*m];
    [λ[i] ~ (α .* (TRM* (β .* (TRM' * E) ./ (TRM' * N))))[i] for i = 1:n];
    [D(S[i]) ~ - λ[i] .* S[i] for i in 1:n];
    [D(E[i]) ~ λ[i] .* S[i] - γₑ * E[i] for i in 1:n];
    [D(R[i]) ~ γₑ * E[i] for i in 1:n]; 
    [D(α[i]) ~ 0. for i in 1:n];
];

@named epi_system = ODESystem(eqs, t)

# Initial conditions 
S0 = [infocomunas[comuna].poblacion for comuna in comunas2];
E0 = 10 * rand(n) .+ 5; # agregué gente random a una comuna random 
R0 = 0.1 * ones(n);
# S0 = [1200., 1550.]; E0 = [3.,5.,]; R0 = [0.,0.]; 

total = S0 + E0 + R0
# Vector de condiciones iniciales 
u0 = [
    [α[i] => 0.1 for i in 1:n];
    [S[i] => S0[i] for i in 1:n]; 
    [E[i] => E0[i] for i in 1:n]; 
    [R[i] => R0[i] for i in 1:n]
];

beta = [1., 5.] # hay que probar qué tan sensible es c/r al segundo valor β₂

p = [
    [γₑ => 1/5.]; # midiendo en semanas... un lío con la matrix P 
    [N[i] => total[i] for i = 1:n];
    [β[i] => beta[i] for i = 1:m];
]
p1 = [
    [γₑ => 1/7]; # midiendo en semanas... un lío con la matrix P 
    [N[i] => total[i] for i = 1:n];
    [β[i] => beta[i] for i = 1:m];
]
# Hay que hacer un collect en algún lado! Chriss Rackauckas

simple_epi_system = structural_simplify(epi_system);


prob = ODEProblem(simple_epi_system, u0, (0.0,400.0), p);
#prob1 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p1);
#prob2 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p);
sol = solve(prob, Tsit5(), saveat = 1.);
#sol1 = solve(prob);
#sol1 =  solve(prob1, Tsit5());

using Plots: plot, plot!

plot(sol, vars = (t, E), legend =:none)
plot(sol, vars = (t, S), legend =:none)
#plot(sol1, vars = (t, E), legend =:none)
#plot(sol, vars = (incidence, t, E, S), legend =:none)

#indexes(i) = ((i-1)*n + 1):(i*n)
indexes(3) # 1 = S, 2 = E, 3 = R 

ts = 0.:1.:300.

#solution(t, i) = sol(t)[indexes(i)]

#plot(sol'[:, indexes(2)]./sol'[1, indexes(1)])

rango = 1:2

a_plot = plot(title = "Incidencias") 
[plot!(a_plot, sol[E[i], :] ./ sol[S[i], 1], label = infocomunas[comunas2[i]].nombre) for i in rango]
display(a_plot)

#=
Jacobiano
=#


jac = generate_jacobian(simple_epi_system); # esto debería ser la función que devuelve el jacobiano 

system_jacobian = eval(jac[1]);

pvec = ModelingToolkit.varmap_to_vars(p,parameters(simple_epi_system));
u0vec = ModelingToolkit.varmap_to_vars(u0,states(simple_epi_system));


system_jacobian(u0vec, pvec, 0.)


rhs_func = eval(generate_function(simple_epi_system)[1]) ; 
rhs_func(u0vec, pvec, 0.)
#incidence(Ei, i::Int) = Ei/S0[i]
#incidence(t, E, S) = (t, E/S)
# Quiero plotear la incidencia... E/S0 

#=
@parameters pₑ pᵢ pₘ γₑ γᵢ γₘ φ 
@variables t S[1:n](t) E[1:n](t) R[1:n](t) λ[1:n](t)
@variables α(t)
@parameters β[1:m]
A = [
    0.8 0.2; 
    0.7 0.3;  
    0.5 0.5; 
    0.6 0.4;
]

TRM = [
    0.8 0.2; 
    0.7 0.3; 
    0.5 0.5; 
    0.6 0.4;
]
 

D = Differential(t) # define an operator for the differentiation w.r.t. time

# your first ODE, consisting of a single equation, indicated by ~
@named fol_separate = ODESystem([ RHS  ~ (1 - x)/τ,
                                  D(x) ~ RHS ])









eqs = [
    [λ[i] ~ control(t) * (TRM* (β .* (TRM' * E) ./ (TRM' * (S + E + R))))[i] for i = 1:n];
    [D(S[i]) ~ - λ[i]* S[i] for i in 1:n];
    [D(E[i]) ~ λ[i]* S[i] - γₑ * E[i] for i in 1:n]
    [D(R[i]) ~ γₑ * E[i] for i in 1:n]
]


control(3)
 #[E => E0], [R => R0])



prob = ODEProblem(structural_simplify(fol_separate), [x => 0.0], (0.0,10.0), [τ => 3.0])
=#
