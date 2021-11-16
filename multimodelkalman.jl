#=
MultipleModelKalmanFilter 
=#

using KalmanFilter 
using ComponentArrays
using LinearAlgebra: I, Diagonal 

# dinámica y jacobiano de la dinámica ya están definidos 
#epi_dynamics(x, α, p, t) 
#epi_jacobian(x, α,  p, t)

generalrk4 = KalmanFilter.GeneralRK4(epi_dynamics, epi_jacobian, dt)

# Ahora un CommonUpdater 
makeF(filter_p) = Diagonal(ModelingToolkit.varmap_to_vars(make_max_vals(filter_p.Fmaxval), states(simple_episys_uknown)))
dims = length(states(simple_episys_uknown))
updater = KalmanFilter.CommonUpdater(dims, generalrk4, makeF, x -> x)


# Un CommonObserver 
obs_exp = [[C[i] for i in 1:n]..., 
        [S[i] + E[i] + I[i] + R[i] for i in 1:n]...
    ]
H = get_observacion_matrix(obs_exp, simple_episys_uknown)
H = SparseMatrixCSC(H)

makeG(filter_p) = Diagonal([filter_p.Gmaxval * ones(n); 10. * ones(n)])
observer = KalmanFilter.CommonObserver(H, makeG, x -> x, 1.)

# Parameters 

priors = [0.1, 0.6, 0.3] 

param = ComponentArray( p = p_realvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [5e-2, 8e-8, 8e-8, 8e-8, 5e-3, 2e-4],
                                    Gmaxval = 100., 
                                    initialcov = [0.1, 0.007, 0.005, 0.07, 0.07, 2e-1]
                                    )
                        )

param2 = ComponentArray( p = p_realvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [1e-2, 2e-9, 2e-8, 2e-8, 1e-3, 8e-5],
                                    Gmaxval = 100., 
                                    initialcov = [0.1, 0.007, 0.005, 0.07, 0.07, 2e-1]
                                    )
                        )

param3 = ComponentArray( p = p_realvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [1e-2, 2e-9, 2e-8, 2e-8, 1e-3, 8e-5],
                                    Gmaxval = 100., 
                                    initialcov = [0.08, 0.0007, 0.0005, 0.007, 0.007, 8e-2]
                                    )
                        )
param4 = ComponentArray( p = pvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [5e-2, 8e-8, 8e-8, 8e-8, 5e-3, 2e-4],
                                    Gmaxval = 100., 
                                    initialcov = [0.1, 0.007, 0.005, 0.07, 0.07, 2e-1]
                                    )
                        )

param5 = ComponentArray( p = pvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [1e-2, 2e-9, 2e-8, 2e-8, 1e-3, 8e-5],
                                    Gmaxval = 100., 
                                    initialcov = [0.1, 0.007, 0.005, 0.07, 0.07, 2e-1]
                                    )
                        )


param6 = ComponentArray( p = p_realvec,
                        filter_p = ComponentArray(
                                    Fmaxval = [1e-2, 2e-9, 2e-8, 2e-8, 1e-3, 8e-5],
                                    Gmaxval = 100., 
                                    initialcov = [0.08, 0.0007, 0.0005, 0.007, 0.007, 8e-2]
                                    )
                        )
param1 = ComponentArray(p = ComponentArray(wn = sqrt(4.)), filter_p = ComponentArray(wn = sqrt(4.)))
param2 = ComponentArray(p = ComponentArray(wn = sqrt(4.4)), filter_p = ComponentArray(wn = sqrt(4.4)))
param3 = ComponentArray(p = ComponentArray(wn = sqrt(4.8)), filter_p = ComponentArray(wn = sqrt(4.8)))

u0 
P0 = Diagonal(ModelingToolkit.varmap_to_vars(make_max_vals(filter_p.initialcov), states(simple_episys_uknown)))

X0 = ComponentArray(x = u0, P = P0)

models = copy([
            KalmanFilter.SimpleKalmanEstimation(param1, copy(X0), copy(X0)), 
            KalmanFilter.SimpleKalmanEstimation(param2, copy(X0), copy(X0)), 
            KalmanFilter.SimpleKalmanEstimation(param3, copy(X0), copy(X0))
        ])

