#= 
Leer movilidad comunal 
=# 

using CSV
using DataFrames
using TimeSeries 


mobfile = "..\\Datos-COVID19-MINSAL\\output\\producto82\\ISCI_weeks.csv" 

dfsalidas = DataFrame(CSV.File(mobfile))

dfsalidasRM = filter(row -> row.region == 13, dfsalidas)

comunas = unique(dfsalidasRM.comuna)
get_salidas(comuna) = filter(row -> row.comuna == comuna, dfsalidasRM).var_salidas  

datamap = Dict(
    comuna => get_salidas(comuna) for comuna in comunas  
) 

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

comunas = collect(keys(datamap))
initial_mob = 0.5 * ones(length(comunas));# get data para esto 

# estoy pensando que esto debería ir al revés... tiempo pal lado 
mob_in_time = Array{Float64, 2}(undef, length(comunas), length(datamap[comunas[1]]))
for (i,comuna) in enumerate(comunas)
    mob_in_time[i, :] = datamap[comuna]
end 
mob_in_time 

comunas2 = comunas[1:2]
mob_in_time2 = Array{Float64, 2}(undef, length(comunas2), length(datamap[comunas[1]]))
for (i,comuna) in enumerate(comunas2)
    mob_in_time2[i, :] = datamap[comuna]
end 
mob_in_time2 
Pt2 = makePmatrix([0.5,0.5], mob_in_time2);


Pt = makePmatrix(initial_mob, mob_in_time);
calculate_semana(t) = floor(Int, t/7) + 1
#index(t) = floor(Int, t) + 1

copyPt = copy(Pt); 
struct DataMatrix
    """
    3d-matrix of dimensions (ncomunas, 2, nsemanas)."""
    Pt::Array{Float64, 3}
end 

function (dm::DataMatrix)(t::Float64, i, j)
    dm.Pt[i,j,calculate_semana(t)]
end 


#dm = DataMatrix(Pt)
dm2 = DataMatrix(Pt2)  
residence_times_matrix2(t, i, j) = dm2(t, i, j)
#@time dm(5.); # 0.000008 seconds (1 allocation: 64 bytes)
#residence_times_matrix(t) = dm(t)
#@time residence_times_matrix(5.); # 0.000022 seconds (4 allocations: 976 bytes)

#=
@code_warntype  residence_times_matrix(5.)
@code_warntype makePmatrix(initial_mob, mob_in_time)

@trace dm(5.); =# 
