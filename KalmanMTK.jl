using KalmanFilter
using LinearAlgebra: Diagonal
using StaticArrays
using SparseArrays

#= 
Dinámica del sistema 

Restricciones y todo eso 
https://academic.csuohio.edu/simond/pubs/IETKalman.pdf
=# 
epi_dynamics(x, α, p, t) = rhs_func(x, p, t)  
epi_jacobian(x, α,  p, t) = system_jacobian(x, p, t) 

#=
Condiciones iniciales 
=# 

dispersion(ν, x) = Diagonal(ν .* x)
a₀ = 0.0135
F = (x) -> dispersion(0.1 * ones(length(x)), x)

max_values = [
    make_alpha(simple_episys_uknown, 5e-3);
    [S[i] => S0[i] for i in 1:n]; 
    [E[i] => 0.07 * S0[i] for i in 1:n];
    [I[i] => 0.07 * S0[i] for i in 1:n]; 
    [R[i] => 0.07 * S0[i] for i in 1:n];
    [C[i] => 0.07 * S0[i] for i in 1:n];
]; 

max_values_vec = ModelingToolkit.varmap_to_vars(max_values, states(simple_episys_uknown))

#=
(S + E + R)./ total


function state_integrity(x, total) 

    x[1:5] = x[1:5] * total/sum(x[1:5])
    x = max.(x, 0.)
    #x[7] = min(x[7], 5.)
    x
end
=#
#= 
Eventualmente le agregaré el índice de α también 
=#
struct StateIntegrity{A <: AbstractArray}
    indexS::A
    indexE::A
    indexR::A 
end 

function (SI::StateIntegrity)(x, total)
    #oldS = x[SI.indexS]; oldE = x[SI.indexE]; oldR = x[SI.indexR]; 
    #suma = oldS + oldE + oldR 
    #x[SI.indexS] .= oldS .* total ./ suma  
    #x[SI.indexE] .= oldE .* total ./ suma  
    #x[SI.indexR] .= oldR .* total ./ suma  
    x .= max.(x, 0.)
end 
# a esto le hace falta testearlo....


# MKT guarda las variables en un orden distinto... 
#https://mtk.sciml.ai/dev/basics/FAQ/#Frequently-Asked-Questions
#indexof(sym,syms) = findfirst(isequal(sym),syms)

getstateindexs(sym_state, system) = [findfirst(isequal(sym_state[i]), states(system)) for i in 1:n] 
indexS = getstateindexs(S, simple_episys_uknown)
indexE = getstateindexs(E, simple_episys_uknown)
indexR = getstateindexs(R, simple_episys_uknown)

#SI1 = StateIntegrity(indexS, indexE, indexR)
SI2 = StateIntegrity(SVector{n}(indexS), SVector{n}(indexE), SVector{n}(indexR))

#@time SI1(u0vec, total);
#@time SI2(u0vec, total);

observation_integrity(x) = max.(x,0.)



lowpass_parameters = ModelingToolkit.varmap_to_vars(
    [
        [S[i] => 1. for i in 1:n]; 
        [E[i] => 1. for i in 1:n]; 
        [R[i] => 1. for i in 1:n]; 
        [I[i] => 1. for i in 1:n]; 
        [C[i] => 1. for i in 1:n]; 
        #[α[i] => 1. for i in 1:n];
        make_alpha(simple_episys_uknown, lowpass_alpha)
        ]
    ,states(simple_episys_uknown)
);

#observaciones = synthetic_obs;
#using Statistics: 

#moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(lenght(vs)-(n-1))]
moving_average(vs,n) = [sum((@view vs[i:(i+n-1), j]))/n for i in 1:(size(vs)[1]-(n-1)), j in 1:size(vs)[2]]


#=begin 
    rkx = KalmanFilter.RK4Dx(epi_dynamics, epi_jacobian, pvec, dt)

    Q = Diagonal(sqrt(dt) * ones(length(u0vec)))

    nlupdater = NLUpdater(rkx, x -> F(max_values_vec), Q, copy(u0vec), 0., 0., x -> SI2(x, total))

    #epi_dynamics(u0vec,0., pvec, 0.) 

    #=
    Observaciones 
    =
    interpolated_data = interpolate_data(confirmedmap[comunas[1]])

    observaciones = Array{Float64, 2}(undef, length(interpolated_data), length(comunas2) + n)
    for (i,comuna) in enumerate(comunas2)
        observaciones[:,i] .= interpolate_data(confirmedmap[comuna])
    end 
    # perfect observations restriction to assure S + E + I + R = N 
    # no es lo más eficiente en términos de memoria....
    
    for i in 1:n 
        observaciones[:, i + length(comunas2)] .= total[i]
    end 
    =# 

    observaciones = synthetic_obs;
    n_obs = size(observaciones)[1]
    for i in 1:n 
        observaciones = [observaciones total[i]*ones(n_obs)]
    end
    system = KalmanFilter.Measurements(observaciones, dt)

    obs_exp = [C[1], C[2], 
        [S[i] + E[i] + I[i] + R[i] for i in 1:n]...  
    ]
    H = get_observacion_matrix(obs_exp, simple_episys_uknown)

    G = ones(2 * n)
    G[1:n] *= 1000.; G[n+1:2n] *= 2.

    observer = KalmanFilter.LinearObserver(H, zeros(2n), G, observation_integrity)
    #=
    Iterator y matriz de covarianzas inicial 
    =# 
    Pfunc = (x) -> F(x) * F(x)'
    iterator = KalmanFilter.LinearKalmanIterator(u0vec, Pfunc(max_values_vec), nlupdater, observer, system, dt, lowpass_parameters)
end =# 
#=
Correr todo 
=#
#results, ensamble = KalmanFilter.full_iteration(iterator, dt, Nmediciones, t -> 0., 1)
#@enter KalmanFilter.full_iteration(iterator, dt, Nmediciones, t -> 0., 1)

#============================================
Loss function
He notado que mientras más cercano a 0. es el 
parámetro del lowpass filter para la tasa de 
contagio, más depende la salida de la condición
inicial. Quiero probar dintintos valores iniciales
para α₀.
=============================================#

obs_exp = [[C[i] for i in 1:n]..., 
        [S[i] + E[i] + I[i] + R[i] for i in 1:n]...
    ]
H = get_observacion_matrix(obs_exp, simple_episys_uknown)
H = SparseMatrixCSC(H)
function loss(analysis, observaciones, rango)
    #sum(abs.(results.analysis[rango, 3:4] - observaciones[rango, :]))
    sum(abs.((analysis * H')[rango,:] - observaciones[rango,:]))
end 

function kalman_iteration(u0, p)
    u0vec = ModelingToolkit.varmap_to_vars(u0,states(simple_episys_uknown));
    pvec = ModelingToolkit.varmap_to_vars(p,parameters(simple_episys_uknown));
    rkx = KalmanFilter.RK4Dx(epi_dynamics, epi_jacobian, pvec, dt)

    Q = Diagonal(sqrt(dt) * ones(length(u0vec)))
    nlupdater = NLUpdater(rkx, x -> F(u0vec), Q, copy(u0vec), 0., 0., x -> SI2(x, total))
    system = KalmanFilter.Measurements(observaciones, dt)

    G = ones(2 * n)
    G[1:n] *= 1000.; G[n+1:2n] *= 10.

    observer = KalmanFilter.LinearObserver(H, zeros(2n), G, observation_integrity)

    #=
    Iterator y matriz de covarianzas inicial 
    =# 
    Pfunc = (x) -> F(x) * F(x)'
    iterator = KalmanFilter.LinearKalmanIterator(u0vec, Pfunc(max_values_vec), nlupdater, observer, system, dt, lowpass_parameters) 

    #results, ensamble = KalmanFilter.full_iteration(iterator, dt, Nmediciones, t -> 0., 1) 
    
    results, ensamble, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs = full_iteration_saver(iterator, dt, Nmediciones, t -> 0., 1)
    xs, Ps = rts_smoother(results, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs, Nmediciones)

    results, xs, Ps
end 

function initial_u0(a0)
    u0 = [
        make_alpha(simple_episys_uknown, a0);
        [S[i] => S0[i] for i in 1:n]; 
        [E[i] => E0[i] for i in 1:n]; 
        [R[i] => R0[i] for i in 1:n];
        [I[i] => I0[i] for i in 1:n];
        [C[i] => C0[i] for i in 1:n];
    ];
end 


function create_p(x)
    @inbounds a0, beta2, gamma_e, gamma_i = x
    p = [
            γₑ => gamma_e, 
            γᵢ => gamma_i, 
            [N[i] => total[i] for i = 1:n]...,
            β[1] => 1.0, 
            β[2] => beta2
        #[γₑ => 1/5.3, γᵢ => 1/9]; # midiendo en semanas... un lío con la matrix P 
        #[γₑ => 1/5.1, γᵢ => 1/8]; # midiendo en semanas... un lío con la matrix P 
    ]
    p
end 


function loss_from_alpha0(a0, p, rango = 1:100)
    u0 = initial_u0(a0)    
    results, xs, Ps = kalman_iteration(u0, p)
    
    println("")
    #loss(results.analysis, observaciones, rango)
    loss(xs, observaciones, rango)
end 


#loss_from_alpha0(a₀, p_real, 1:400)
#loss_from_alpha0(a₀, create_p([0.016728510433320562, 34.87213978597166,  0.1451405896110646,  0.16229284534863325]), 1:400)


#@enter kalman_iteration(initial_u0(a₀), p);
#results, xs, Ps = kalman_iteration(initial_u0(best_candidate(opt)[1]), create_p(best_candidate(opt)));
#run kalman_iteration(initial_u0(a₀), p)

#@run kalman_iteration(initial_u0(a₀), p)

#= 
a0s = 0.001:0.001:0.01
a0s = 1e-7:1e-7:1e-6

loss_from_alpha0(1e-10) 

losses_a0s = loss_from_alpha0.(a0s)  
=#
#using Plots: scatter, scatter!
#scatter(a0s, losses_a0s)

#============================================
Plotting 
=============================================#

# Agregar un campo opcional con nombres
# (por defecto es el nombre de las mismas variables)


#=plot_smoothed!(a_plot, ts, xs, 1., i, subplot = state) # 10^5 
#plot_smoothed!(a_plot, ts, xs, 1., 2, 1)
#plot_scnotation!(a_plot, sol, S, 1)

plot_smoothed!(a_plot, ts, xs, 1., n + i, ) # 10^2
#plot_smoothed!(a_plot, ts, xs, 1., 4, 2)
#plot_scnotation!(a_plot, sol, E, 2)

plot_smoothed!(a_plot, ts, xs, 1., 2n + i, 3) # 10^2
#plot_smoothed!(a_plot, ts, xs, 1., 6, 3)
#plot_scnotation!(a_plot, sol, I, 3)

plot_smoothed!(a_plot, ts, xs, 1., 3n + i, 4) # 10^4
#plot_smoothed!(a_plot, ts, xs, 1., 8, 4)
#plot_scnotation!(a_plot, sol, R, 4)

plot_smoothed!(a_plot, ts, xs, 1., 4n + i, 5) # 10{4}
#plot_smoothed!(a_plot, ts, xs, 1., 10, 5)
#plot_scnotation!(a_plot, sol, C, 5)

plot_smoothed!(a_plot, ts, xs, 1., 5n + i, 6) # 1
#plot_smoothed!(a_plot, ts, xs, 1., 12, 6) # 1 =#

#plot!(ts, control_pieces.(ts), ylabel = "α(t)", subplot = 6, legend=:none, xlabel = "t")
#savefig(folder * "kalmana_llstates" * make_img_name(p_real) * ".svg")


#plot!(ts, xs[:,i]./total[i], ribbon = sqrt.(Ps[i,i,:])./total[i])

#=
plot(results, ts, 1)
plot(results, ts, 1, ylims = (0.,1.9e5), title = "Susceptibles")

plot!(sol, vars = (t,S))

plot(results, ts, 3, title = "Expuestos")
plot!(results, ts, 4)

plot!(ts, xs[:,4], label = "s2")
plot!(sol, vars = (t,E))

plot(results, ts, 5, title = "Infectados")
plot!(results, ts, 6)
plot!(sol, vars = (t,I))

plot(results, ts, 7, title = "Recuperados")
plot!(results, ts, 8) 
plot!(sol, vars = (t,R))

plot(results, ts, 9, title = "Acumulados")
plot!(results, ts, 10)
plot!(sol, vars = (t,C))
plot!(ts, xs[:,9], label = "sC1")
plot!(ts, xs[:,10], label = "sC2")
plot!(ts, observaciones[rango, 1], label = "C[1] obs") =#
#plot!(ts, observaciones[rango, 2], label = "C[2] obs") 

#plot(results, ts, 11, title = "Tasa de contagios")
#plot!(results, ts,12) 
#plot!(ts, xs[:,11], label = "sα1")
#plot!(ts, xs[:,12], label = "sα2")
#plot!(ts, control_pieces.(ts))

#=
find_in_states(sym, system) = findfirst(isequal(sym),states(system))

function result_dic(analysis, system) 
    results_from_sym(sym) = analysis[:,find_in_states(sym, system)]
    dict = [state => results_from_sym(state) for state in states(system)]
end 

function vec_from_dic(dic, system)
    hcat(ModelingToolkit.varmap_to_vars(dic, states(system))...)'
end 
=#

#=
"""
- `sym_expression`: using variables from states(system). 
- `ordered_results`: an 2d array 
        filas: states, in the order given by ModelingToolkit.varmap_to_vars.
        columns: timestamps 
# Example 
```julia
julia> eval_expression(C[1] + E[2] + 5R[1], ordered_results, simple_episys_uknown)
```
"""
function eval_expression(sym_expression, ordered_results, system)
    number_of_states, total_steps = size(ordered_results)
    func = eval(build_function(sym_expression, states(system),  expression=Val{false}))
    
    evaluated_array = Vector{Float64}(undef, total_steps)
    for k in 1:total_steps
        states_at_k = ordered_results[:,k]
        evaluated_array[k] = func(states_at_k)
    end 
    evaluated_array
end 

function evaluate_expression_in_analysis(exp, analysis, system)
    dict = result_dic(analysis, system)
    eval_expression(exp, vec_from_dic(dict, system), system)
end 

total_class_1 = evaluate_expression_in_analysis(S[1] + E[1] + R[1] + I[1], results.analysis, simple_episys_uknown);
total_class_1_rts = evaluate_expression_in_analysis(S[1] + E[1] + R[1] + I[1], xs, simple_episys_uknown);

plot(ts, total_class_1, ylims = (136300.,136450.))
plot!(ts, total_class_1_rts, label = "rts")
=#
#=
plot(results, ts, 3, title = "Expuestos")
plot!(ts, xs[:,3])

plot(results, ts, 1, title = "Susceptibles")
plot!(ts, xs[:,1])

plot(results, ts, 5, title = "Infected")
plot!(ts, xs[:,5])

plot(results, ts, 7, title = "Recovered")
plot!(ts, xs[:,7])

plot(results, ts, 9, title = "Cumulated")
plot!(ts, xs[:,9])


plot(results, ts, 11, title = "Control")
plot!(ts, xs[:,11])

=#

#=
if i == 6
        @series begin
            seriestype := :path
            label --> "Observaciones"
            rango_ts = (rango[1]+1):rango[end]
            rango_obs = rango[1]:(rango[end]-1)
            ts[rango_ts], r.observations[rango_obs]
        end
    end
=#

#=
En el caso sintético se obtienen buenos resultados!
Creo que necesito corregir los casos acumulados... están dando problema...
supongo que ese salto raro al cambiar las mediciones. 
Si logro hacer optimización de parámetros... podría funcionar, es claro que 
estos parámetros no son los mejores, el forzante creo que apaña a eso. 
=#