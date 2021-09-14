using KalmanFilter
using LinearAlgebra: Diagonal
using StaticArrays

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

u0vec
dispersion(ν, x) = Diagonal(ν .* x)

max_values = [
    [α[i] => 1e-2 for i in 1:n];
    [S[i] => S0[i] for i in 1:n]; 
    [E[i] => 0.3 * S0[i] for i in 1:n];
    [I[i] => 0.3 * S0[i] for i in 1:n]; 
    [R[i] => 0.3 * S0[i] for i in 1:n];
    [C[i] => 0.3 * S0[i] for i in 1:n];
]; 

max_values_vec = ModelingToolkit.varmap_to_vars(max_values, states(simple_epi_system))

dt = 1.#0.1 # Intervalo de sampleo de observaciones 
T = 400. # Se resolverá el problema en el intervalo [0,T]
Nmediciones = Int(T/dt) # número de mediciones 
ts = 0.:dt:(T-dt) # grilla de tiempos 
 
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

getstateindexs(sym_state) = [findfirst(isequal(sym_state[i]), states(epi_system)) for i in 1:n] 
indexS = getstateindexs(S)
indexE = getstateindexs(E)
indexR = getstateindexs(R)

#SI1 = StateIntegrity(indexS, indexE, indexR)
SI2 = StateIntegrity(SVector{n}(indexS), SVector{n}(indexE), SVector{n}(indexR))

#@time SI1(u0vec, total);
@time SI2(u0vec, total);

observation_integrity(x) = max.(x,0.)


dt = 1. 

lowpass_parameters = ModelingToolkit.varmap_to_vars(
    [
        [S[i] => 0.8 for i in 1:n]; 
        [E[i] => 0.8 for i in 1:n]; 
        [R[i] => 0.8 for i in 1:n]; 
        [I[i] => 0.8 for i in 1:n]; 
        [C[i] => 0.8 for i in 1:n]; 
        [α[i] => 0.3 for i in 1:n];
        ]
    ,states(simple_epi_system)
);

begin 
    rkx = KalmanFilter.RK4Dx(epi_dynamics, epi_jacobian, pvec, dt)

    Q = Diagonal(sqrt(dt) * ones(length(u0vec)))

    nlupdater = NLUpdater(rkx, x -> F(max_values_vec), Q, copy(u0vec), 0., 0., x -> SI2(x, total))

    #epi_dynamics(u0vec,0., pvec, 0.) 

    #=
    Observaciones 
    =#
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
    system = KalmanFilter.Measurements(observaciones, dt)

    obs_exp = [C[1], C[2], 
        [S[i] + E[i] + I[i] + R[i] for i in 1:n]...  
    ]
    H = get_observacion_matrix(obs_exp, simple_epi_system)

    G = ones(2 * n)
    G[1:n] *= 1000.; G[n+1:2n] *= 2.

    observer = KalmanFilter.LinearObserver(H, zeros(2n), G, observation_integrity)
    #=
    Iterator y matriz de covarianzas inicial 
    =# 
    Pfunc = (x) -> F(x) * F(x)'
    iterator = KalmanFilter.LinearKalmanIterator(u0vec, Pfunc(max_values_vec), nlupdater, observer, system, dt, lowpass_parameters)
end 
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

function loss(results, observaciones, rango)
    sum(abs.(results.analysis[rango, 3:4] - observaciones[rango, :]))
end 

function kalman_iteration(u0, p)
    u0vec = ModelingToolkit.varmap_to_vars(u0,states(simple_epi_system));
    pvec = ModelingToolkit.varmap_to_vars(p,parameters(simple_epi_system));
    rkx = KalmanFilter.RK4Dx(epi_dynamics, epi_jacobian, pvec, dt)

    Q = Diagonal(sqrt(dt) * ones(length(u0vec)))
    nlupdater = NLUpdater(rkx, x -> F(u0vec), Q, copy(u0vec), 0., 0., x -> SI2(x, total))
    system = KalmanFilter.Measurements(observaciones, dt)


    obs_exp = [C[1], C[2]]
    H = get_observacion_matrix(obs_exp, simple_epi_system)

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
        [α[i] => a0 for i in 1:n];
        [S[i] => S0[i] for i in 1:n]; 
        [E[i] => E0[i] for i in 1:n]; 
        [R[i] => R0[i] for i in 1:n];
        [I[i] => I0[i] for i in 1:n];
        [C[i] => C0[i] for i in 1:n];
    ];
end 

function loss_from_alpha0(a0, rango = 1:100)
    u0 = initial_u0(a0)    
    results = kalman_iteration(u0, p)
    
    println("")
    loss(results, observaciones, rango)
end 


results, xs, Ps = kalman_iteration(initial_u0(a₀), p)
#@run kalman_iteration(initial_u0(a₀), p)

#= 
a0s = 0.001:0.001:0.01
a0s = 1e-7:1e-7:1e-6

loss_from_alpha0(1e-10) 

losses_a0s = loss_from_alpha0.(a0s)  
plot(a0s, losses_a0s)

#============================================
Plotting 
=============================================#

rango = 1:length(ts)

using Plots: plot, plot!
find_in_states(sym, system) = findfirst(isequal(sym),states(system))

function result_dic(analysis, system) 
    results_from_sym(sym) = analysis[:,find_in_states(sym, system)]
    dict = [state => results_from_sym(state) for state in states(system)]
end 

dict = result_dic(results.analysis, simple_epi_system)

function vec_from_dic(dic, system)
    hcat(ModelingToolkit.varmap_to_vars(dic, states(system))...)'
end 


"""
- `sym_expression`: using variables from states(system). 
- `ordered_results`: an 2d array 
        filas: states, in the order given by ModelingToolkit.varmap_to_vars.
        columns: timestamps 
# Example 
```julia
julia> eval_expression(C[1] + E[2] + 5R[1], ordered_results, simple_epi_system)
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

total_class_1 = eval_expression(S[1] + E[1] + R[1] + I[1], vec_from_dic(dict, simple_epi_system), simple_epi_system)

plot(ts, total_class_1)