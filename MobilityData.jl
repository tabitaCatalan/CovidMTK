#= 
Leer movilidad comunal 
=# 

using CSV
using DataFrames
using TimeSeries 

#============================================================================
Definición de funciones para procesar los datos 
============================================================================# 

"""
    process_mob_csv(file)
Devuelve un Vector con los códigos de las comunas de Santiago y un diccionario asociando 
cada código a su respectiva movilidad.
- `file`: path al archivo de datos de mobilidad (producto 82)
```
julia> comunas
51-element Vector{Int64}:
 13111
     ⋮
 13112
 julia> datamap
 Dict{Int64, Vector{Float64}} with 51 entries:
   13111 => [0.992765, 1.00055,  …  0.871758, 0.877422, 0.897… 
   13503 => [0.990494, 0.992677,  …  0.892761, 1.00352, 0.9580…  
   ⋮     => ⋮
```
"""
function process_mob_csv(file)
    dfsalidas = DataFrame(CSV.File(mobfile))

    dfsalidasRM = filter(row -> row.region == 13, dfsalidas)

    comunas = unique(dfsalidasRM.comuna)
    get_salidas(comuna) = filter(row -> row.comuna == comuna, dfsalidasRM).var_salidas  

    datamap = Dict(
        comuna => get_salidas(comuna) for comuna in comunas  
    ) 
    comunas, datamap
end 


make_initial_homogeneous_mob(initial_frac_home_time, number_of_comunas) = initial_frac_home_time * ones(number_of_comunas)

"""
    make_mob_in_time(comunas)
Construye una matriz de mobilidad semanal para todas las comunas elegidas 
- `comunas::Array{Int}` lista de códigos de las comunas 
"""
function make_mob_in_time(selected_comunas, datamap)
    mob_in_time = Array{Float64, 2}(undef, length(selected_comunas), length(datamap[selected_comunas[1]]))
    for (i,comuna) in enumerate(selected_comunas)
        mob_in_time[i, :] .= datamap[comuna]
    end 
    mob_in_time 
end 

"""
Devuelve una matriz de tiempos de residencia que depende del tiempo 
- `initial_mob`: array_(i), donde i es la comuna 
- `mob_in_time`: array_(i,j), donde i es la comuna y j la semana, variaciones
    con respecto a la cantidad inicial 
"""
function makePmatrix(initial_mob, mob_in_time)
    ncomunas, nsemanas = size(mob_in_time)
    Pt = Array{Float64, 3}(undef, ncomunas, 2, nsemanas) 
    for semana in 1:nsemanas 
        variation_t = @view mob_in_time[:, semana]
        Pt[:,:,semana] = [(1 .- variation_t .* initial_mob) (variation_t .* initial_mob)]
    end 
    Pt  
end  

calculate_semana(t) = floor(Int, t/7) + 1

struct DataMatrix
    """
    3d-matrix of dimensions (ncomunas, 2, nsemanas)."""
    Pt::Array{Float64, 3}
end 

function (dm::DataMatrix)(t::Float64, i, j)
    dm.Pt[i,j,calculate_semana(t)]
end 


#============================================================================
Ejecutar las funciones 
============================================================================# 

if Sys.islinux()
    mobfile = "../Datos-COVID19-MINSAL/output/producto82/ISCI_weeks.csv" 
else 
    mobfile = "..\\Datos-COVID19-MINSAL\\output\\producto82\\ISCI_weeks.csv"
end

comunas, datamap = process_mob_csv(mobfile) 

try
    # Si es que initial_data_mob está definido, usarlo. S
    global initial_mob = initial_data_mob 
catch # si no está definido, usar mobilidad inicial homogenea
    global initial_mob = make_initial_homogeneous_mob(initial_frac_home_time, length(comunas2))
end 


#@time dm(5.); # 0.000008 seconds (1 allocation: 64 bytes)
#residence_times_matrix(t) = dm(t)
#@time residence_times_matrix(5.); # 0.000022 seconds (4 allocations: 976 bytes)

#= código para ubicar dónde se están haciendo allocations, qué hace que la función demore 
@code_warntype  residence_times_matrix(5.)
@code_warntype makePmatrix(initial_mob, mob_in_time)

@trace dm(5.); =# 
