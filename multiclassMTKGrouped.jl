using ModelingToolkit
using OrdinaryDiffEq
using Plots: plot, plot!

# GROUPED COMUNAS 
n = length(totales_por_clase) # number of social classes 
m = 2 # number of environments 

# register data functions 
@variables t 
@register control_pieces(t)
@register residence_times_matrix(t, i, j)
@parameters N[1:n]
@variables S[1:n](t) E[1:n](t) I[1:n](t) R[1:n](t) C[1:n](t) λ[1:n](t)


if variable_rate 
    @variables γₑ(t) γᵢ(t)
else 
    @parameters γₑ γᵢ
end 

β = get_beta(variable_beta, t) 
#@variables TRM[1:n, 1:m](t)
D = Differential(t) 
#= +28 ... para ajustar los tiempos
mobility data starts 2020-03-02 (epi-week 9)
and confirmed data starts on 2020-03-30 
confirmed prod 15 starts 2020-02-22 
#semana_to_lastday(1) + Day(37)
=#
adjusted_rtm(t,i,j) = residence_times_matrix(t+28, i, j)

# epi_model_unknown_input y epi_model_known_input están definidas en multiclass_model.jl 
@named episys_uknown = epi_model_unknown_input(t, n, m, adjusted_rtm, false, variable_rate, gamma_e_real, gamma_i_real, variable_beta, beta_exterior_real)
@named episys_known = epi_model_known_input(t, n, m, adjusted_rtm, control_pieces, false, variable_rate,gamma_e_real, gamma_i_real, variable_beta, beta_exterior_real)

#=
Initial conditions 
=#

S0 = totales_por_clase
#S0 = [infocomunas[comuna].poblacion for comuna in comunas2];
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

if variable_rate
    global u0 = [u0; make_rate(gamma_e, gamma_i)]
end 
if variable_beta
    global u0 = [u0; make_beta(beta_exterior)]
end 


beta_real = [1., 50.]
beta = [1., beta_exterior] # hay que probar qué tan sensible es c/r al segundo valor β₂


p_real = make_p(episys_known, 1/5.3, 1/8.3, total, beta_real)



p = make_p(episys_uknown, gamma_e, gamma_i, total, beta) # otros valores 

#=p1 = [
    [γₑ => 1/7]; # midiendo en semanas... un lío con la matrix P 
    [N[i] => total[i] for i = 1:n];
    [β[i] => beta[i] for i = 1:m];
]=#

# Synthetic data
#prob_known = ODEProblem(simple_episys_known, u0_real, (0.0,T), p_real, jac = true, sparse = true);

#prob1 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p1);
#prob2 = ODEProblem(simple_epi_system, u0, (0.0,400.0), p);
#sol = solve(prob_known, Tsit5(), saveat = 1.);
#synthetic_obs = [sol[C[1]] sol[C[2]]];


using Plots: plot, plot!, savefig

#=
plot(ts, control_pieces.(ts), )
savefig(folder * "control" * make_img_name(p_real) * ".svg")

plot(sol, vars = (t, E), title = "Exposed", ylabel = "people")
savefig(folder * "exposed" * make_img_name(p_real) * ".svg")
plot!(sol, vars = (t, I), title = "Infected and Exposed", ylabel = "people")
savefig(folder * "infected" * make_img_name(p_real) * ".svg")
plot(sol, vars = (t, S), title = "Susceptibles", ylabel = "people")#, ylim = (0,3.4e5))
savefig(folder * "susceptibles" * make_img_name(p_real) * ".svg")
plot(sol, vars = (t, R), title = "Recovered", ylabel = "people")
savefig(folder * "recovered" * make_img_name(p_real) * ".svg")

plot!(sol, vars = (t, C), title = "Cumulated", ylabel = "people")
savefig(folder * "cumulated" * make_img_name(p_real) * ".svg")
=#
#=
a_plot = plot(sol, vars = (t, C))
plot!(a_plot, ts, observaciones[1:length(ts), 1], label = "data C[1]")
plot!(a_plot, ts, observaciones[1:length(ts), 2], label = "data C[2]")
display(a_plot)

plot(ts, observaciones[1:length(ts), 1], label = "data C[1]", yscale = :log10)
plot(sol, vars = (t, C), title = "Cumulated cases")

plot(ts, observaciones[1:length(ts), 1], label = "data C[1]")
plot!(ts, observaciones[1:length(ts), 2], label = "data C[2]")

#plot(sol1, vars = (t, E), legend =:none)
#plot(sol, vars = (incidence, t, E, S), legend =:none)
=#
#indexes(i) = ((i-1)*n + 1):(i*n)

#solution(t, i) = sol(t)[indexes(i)]

#plot(sol'[:, indexes(2)]./sol'[1, indexes(1)])

rango = 1:2


#=
Better plot, shared x-axis 

a_plot = plot(layout=(4,1),framestyle=:box, link = :x, size = (400, 600))
plot_scnotation!(a_plot, sol, S, 1)
plot_scnotation!(a_plot, sol, E, 2)
plot_scnotation!(a_plot, sol, I, 2)
plot_scnotation!(a_plot, sol, R, 3)
plot_scnotation!(a_plot, sol, C, 3)
plot!(ts, control_pieces.(ts), ylabel = "α(t)", subplot = 4, legend=:none, xlabel = "t")
savefig(folder * "allstates" * make_img_name(p_real) * ".svg")

## mobility 
plot(100 * datamap[comunas2[1]], ylabel = "% variation w/r to initial mobility", label = "municipality 1")
plot!(100 * datamap[comunas2[2]], label = "municipality 2", ylims = (0,100), legend = :bottomright)
savefig(folder * "mobility_$(comunas2[1])_$(comunas2[2]).svg")

a_plot = plot(title = "Incidencias") 
[plot!(a_plot, sol[E[i], :] ./ sol[S[i], 1], label = "municipality $i") for i in rango]
display(a_plot)
savefig(a_plot, folder * "incidence" * make_img_name(p_real) * ".svg")
=# 
#=
Jacobian to use with Kalman Filter 
=#


jac2 = generate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
#jac = generate_jacobian(simple_episys_uknown); # esto debería ser la función que devuelve el jacobiano 

#jac2 = generate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
system_jacobian = eval(jac2[1]);

#jac2 = calculate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
#jac1 = calculate_jacobian(simple_episys_uknown);
pvec = ModelingToolkit.varmap_to_vars(p,parameters(simple_episys_uknown));
u0vec = ModelingToolkit.varmap_to_vars(u0,states(simple_episys_uknown));
#class_names_vec = ModelingToolkit.varmap_to_vars(class_names,states(simple_episys_uknown));

#system_jacobian(u0vec, pvec, 0.)
#build_function(jac1, states(simple_episys_known));

rhs_func = eval(generate_function(simple_episys_uknown)[1]) ; 
rhs_func(u0vec, pvec, 0.)
#incidence(Ei, i::Int) = Ei/S0[i]
#incidence(t, E, S) = (t, E/S)
# Quiero plotear la incidencia... E/S0 
