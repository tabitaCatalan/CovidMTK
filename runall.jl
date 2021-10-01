#a = 1
include("MobilityData.jl")
include("ConfirmedComunas.jl")
include("InterpolateConfirmed.jl")
include("ConfirmedProd15.jl")
include("common_defs.jl")
include("control.jl")
include("multiclassmodel.jl")
include("plotting.jl")
include("multiclassMTK.jl")
include("linear_coeff.jl")
include("smoother.jl")
include("KalmanMTK.jl") 

#include("optim.jl")]

tsdate = Date(2020,03,30):Day(1):Date(2021,04,1) 

