# Smoother 


function full_iteration_saver(iterator::KalmanFilter.LinearKalmanIterator,
    dt, N, control_function, ensamble_size
    )

    dimensions = KalmanFilter.kalman_size(iterator.observer)[1]
    Pnp1n_matrixs = Array{Float64, 3}(undef, dimensions, dimensions, N)
    Pnn_matrixs = Array{Float64, 3}(undef, dimensions, dimensions, N)
    Fn_matrixs = Array{Float64, 3}(undef, dimensions, dimensions, N)
    results = KalmanFilter.InnerStateSeries(N, dimensions)
    ensamble = KalmanFilter.EnsamblesStoring(ensamble_size, dimensions, N)

    for i in 1:N
        #print(i)
        control = control_function(i * dt)
        #observation = KalmanFilter.observe_inner_system(iterator)
        KalmanFilter.add_ensamble!(ensamble, i, iterator)
        #add_state!(results, i, get_inner_state(iterator.system))
        KalmanFilter.add_analysis!(results, i, KalmanFilter.hatx(iterator))


        Xₙ₊₁ₙ = KalmanFilter.forecast(iterator, control)
        KalmanFilter.add_forecast!(results, i, Xₙ₊₁ₙ.x)
        KalmanFilter.add_error_analysis!(results, i, [KalmanFilter.hatP(iterator)[j,j] for j in 1:dimensions])
        #add_error_forecast!(results, i, [forecastP[j,j] for j in 1:dimensions])
        Fn_matrixs[:,:,i] .= KalmanFilter.Mn(iterator.updater)

        Pnp1n_matrixs[:,:,i] .= KalmanFilter.next_hatP(iterator)
        Pnn_matrixs[:,:,i] .= KalmanFilter.hatP(iterator)      

        KalmanFilter.add_observation!(results, i, KalmanFilter.next_iteration!(iterator, control)[1])
        print(Int(KalmanFilter.check_observability(iterator)))
    end
    results, ensamble, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs
end


function rts_smoother(results, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs, N)
    # initialize
    xs = similar(results.analysis)
    xs[end,:] = results.analysis[end,:] # last x 

    Ps = similar(Pnp1n_matrixs)
    Ps[:,:,end] = Pnn_matrixs[:,:,N]

    for n in N-1:-1:1 
        Pnp = Pnn_matrixs[:,:,n]
        #Gnp1 = Cn_matrixs[:,:,n+1] * inv(Pnp1n_matrixs[:,:,n+1]) # la versión continuo discreta solo cambia la matriz de ganacia 
        Fn = Fn_matrixs[:,:,n]
        Gnp1 = Pnp * Fn' * inv(Pnp1n_matrixs[:,:,n+1])
        Ps[:,:,n] = Pnp - Gnp1 * (Pnp1n_matrixs[:,:,n+1] - Ps[:,:,n+1]) * Gnp1'
        xs[n,:] = results.analysis[n,:] + Gnp1 * (xs[n+1,:] - results.forecast[n+1,:])
    end
    xs, Ps
end

#results, ensamble, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs = full_iteration_saver(iterator, dt, Nmediciones, t -> 0., 1)
#xs, Ps = rts_smoother(results, Pnp1n_matrixs, Pnn_matrixs, Fn_matrixs, Nmediciones)

#=
using Plots: plot, plot!
plot(results, ts, 2, error = false)

plot(ts, results.analysis[:,1])
plot!(ts, xs[:,1], label = "RTS smoother")
#plot!(ts, sdesol'[1:end-1, 2], label = "Sol real")
plot!(ts[end-1:end], xs[end-1:end,1]) 

=#