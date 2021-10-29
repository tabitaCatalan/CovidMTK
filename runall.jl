
# grouping_of_municipalities = 3 # high, medium, low 
number_of_selected_municipalities = 5; # con cuántas municipalidades trabajar 
initial_frac_home_time = 0.5
lag_days_confirmed = 7 # número de días para el promedio móvil del número de confirmados 
one_control = false
lowpass_alpha = 0.99
beta_exterior = 60.
number_of_groups = 5 
# grouped = false # No voy a usar esta variable, separaré los archivos importantes en dos 
#=
Voy a separar el código en 3 secciones: 
- Load data: carga todos los datos disponibles y los prepara para ser elegidos
- Data selection and preparation: depende de parametros como el número de comunas elegidas.
    Selecciona los datos requeridos desde los datos completos y los prepara (interpolando
    por ejemplo) para ser usados.
- Common definitions: dt, T, etc. Control used for generating synthetic data. 
- Model and variables: crea las variables simbólicas y el model en MTK 
=#
include("ConfirmedComunas.jl")
include("sort_group_muni_by_ips.jl")
include("MobilityData.jl")


include("InterpolateConfirmed.jl")
include("ConfirmedProd15.jl")

prod15map, lastepiday = process_prod_15(prod15file) 

#=
Ver tabitaCatalan/CovidMTK#1 
Están del más bajo al más alto, al revés 
julia> initial_res_time_matrix #quitando 7/24 de hogar 
5×2 Matrix{Float64}:
 0.672939  0.327061
 0.666404  0.333596
 0.671306  0.328694
 0.663146  0.336854
 0.621711  0.378289

julia> initial_res_time_matrix_with_sleep
5×2 Matrix{Float64}:
 0.76915   0.23085
 0.766318  0.233682
 0.768625  0.231375
 0.762551  0.237449
 0.732865  0.267135
=#
#inital_data_mob = [0.378289, 0.336854, 0.328694, 0.333596, 0.327061] # invertida, descontando 7/24
initial_data_mob = [0.267135, 0.237449, 0.231375, 0.233682, 0.23085] # invertida
initial_mob = initial_data_mob  # en caso de tener una mobilidad inicial, e.g. dada por la encuesta origen destino 

#initial_mob = make_initial_homogeneous_mob(initial_frac_home_time, number_of_groups)
#initial_mob = make_initial_homogeneous_mob(initial_frac_home_time, number_of_selected_municipalities )


include("grouped.jl") 
#include("not-grouped.jl")

interpolated_prod15_grouped = make_grouped_prod15_confirmed(groups, prod15map, lastepiday)
#interpolated_prod15 = make_prod15_confirmed(comunas2, prod15map, lastepiday)

dm = DataMatrix(Pt)
residence_times_matrix(t, i, j) = dm(t, i, j)

#=---------------------------------=#

include("common_defs.jl")
include("control.jl")
include("multiclassmodel.jl")
include("plotting.jl")


#===incluir solo uno de los dos ===#
include("multiclassMTKGrouped.jl")
#include("multiclassMTK.jl")
#=---------------------------------=#

include("linear_coeff.jl")
include("smoother.jl")
#include("smoother-cd.jl")

#moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(lenght(vs)-(n-1))]
moving_average(vs,n) = [sum((@view vs[i:(i+n-1), j]))/n for i in 1:(size(vs)[1]-(n-1)), j in 1:size(vs)[2]]


#===incluir solo uno de los dos ===#
observaciones = moving_average(interpolated_prod15_grouped, lag_days_confirmed)
n_obs = size(observaciones)[1]
for i in 1:n 
    observaciones = [observaciones totales_por_clase[i]*ones(n_obs)]
end


observaciones = moving_average(interpolated_prod15, lag_days_confirmed)
n_obs = size(observaciones)[1]
for i in 1:n 
    observaciones = [observaciones total[i]*ones(n_obs)]
end
#=---------------------------------=#

include("KalmanMTK.jl") 

#include("optim.jl")

#tsdate = Date(2020,03,30):Day(1):Date(2021,04,1) 

a = 1
#=
using ModelingToolkit 
@variables t x(t) y(t) 
D = Differential(t) 
eqs = [D(x) ~ -y , D(y) ~ 2x ] 
@named sys = ODESystem(eqs, t) 
jac = eval(generate_jacobian(sys)[1]);
jac2 = eval(generate_jacobian(sys, sparse = true)[1]);

jac2([2.,1.], [], 3.)
=#

function buscar_id(nombre_comuna)
    filter(comuna -> last(comuna).nombre == nombre_comuna, infocomunas)    
end  

buscar_id("Puente Alto")


