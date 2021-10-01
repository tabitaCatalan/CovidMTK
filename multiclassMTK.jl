using ModelingToolkit, OrdinaryDiffEq
using Plots: plot, plot!


n = length(comunas2) # number of classes 
m = 2 # number of environments 

# register data functions 
@variables t 
@register control_pieces(t)
@register residence_times_matrix2(t, i, j)
@parameters γₑ γᵢ β[1:m] N[1:n]
@variables S[1:n](t) E[1:n](t) I[1:n](t) R[1:n](t) C[1:n](t) λ[1:n](t)
@variables TRM[1:n, 1:m](t)
D = Differential(t) 
#= +28 ... para ajustar los tiempos
mobility data starts 2020-03-02 (epi-week 9)
and confirmed data starts on 2020-03-30 
confirmed prod 15 starts 2020-02-22 
#semana_to_lastday(1) + Day(37)
=#
adjusted_rtm(t,i,j) = residence_times_matrix2(t+28, i, j)

# epi_model_unknown_input y epi_model_known_input están definidas en multiclass_model.jl 
@named episys_uknown = epi_model_unknown_input(t, n, m, adjusted_rtm, false)
@named episys_known = epi_model_known_input(t, n, m, adjusted_rtm, control_pieces)

#=
Initial conditions 
=#


S0 = [infocomunas[comuna].poblacion for comuna in comunas2];
E0 = 10 * rand(n) .+ 20; # agregué gente random a una comuna random 
I0 = 0.1 * ones(n);
R0 = 0.1 * ones(n);
C0 = 0.1 * ones(n);
# S0 = [1200., 1550.]; E0 = [3.,5.,]; R0 = [0.,0.]; 

total = S0 + E0 + R0

# Structural simplifications to eliminate extra variables 
simple_episys_known = structural_simplify(episys_known); # to use in kalman filter 
simple_episys_uknown = structural_simplify(episys_uknown); # to generate synthetic data 


# 0132, 135 es mucho
#=
#Lo ideal sería que funcionara esto........
=#
u0_real = make_x_k(episys_known, S0, E0, I0, R0, C0) 
#=
u0_real = [
    [S[i] => S0[i] for i in 1:n];
    [E[i] => E0[i] for i in 1:n];
    [I[i] => I0[i] for i in 1:n];
    [R[i] => R0[i] for i in 1:n];
    [C[i] => C0[i] for i in 1:n];
] =#

u0 = make_x_uk(episys_uknown, 0.0110, S0, E0, I0, R0, C0)

beta_real = [1., 50.]
beta = [1., 60.] # hay que probar qué tan sensible es c/r al segundo valor β₂


p_real = make_p(episys_known, 1/5.3, 1/8.3, total, beta_real)


p = make_p(episys_uknown, 1/5.1, 1/7.2, total, beta) # otros valores 

#=p1 = [
    [γₑ => 1/7]; # midiendo en semanas... un lío con la matrix P 
    [N[i] => total[i] for i = 1:n];
    [β[i] => beta[i] for i = 1:m];
]=#

# Synthetic data
prob_known = ODEProblem(simple_episys_known, u0_real, (0.0,T), p_real);

#prob1 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p1);
#prob2 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p);
sol = solve(prob_known, Tsit5(), saveat = 1.);
synthetic_obs = [sol[C[1]] sol[C[2]]];


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

#a_plot = plot(title = "Incidencias") 
#[plot!(a_plot, sol[E[i], :] ./ sol[S[i], 1], label = infocomunas[comunas2[i]].nombre) for i in rango]
#display(a_plot)

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