#=
confirmed using product 15 (por semana epidemiológica)
=#

using CSV 
using TimeSeries
using DataFrames 
if Sys.islinux()
    prod15file = "../Datos-COVID19-MINSAL/output/producto15/FechaInicioSintomas.csv"
else 
    prod15file = "C:\\Users\\Tabita\\Documents\\Covid\\Datos-COVID19-MINSAL\\output\\producto15\\FechaInicioSintomas.csv"
end

"""
Entrega el último día de la semana epidemiológica pedida 
# Ejemplo 
Obtener el último día de la semana epidemiológica 6 (SE6)
```julia 
julia> semana_to_lastday(6)
2020-03-28
```
"""
semana_to_lastday(week) = Date(2020,2,22) + (week - 1) * Day(7)

function process_prod_15(file)
    dfprod15 = DataFrame(CSV.File(file))

    dfprod15RM = filter(row -> ismissing(row["Codigo region"]) ? false : row["Codigo region"] .== 13. , dfprod15)

    prod15map = Dict(
        # código comuna => array de confirmados de esa comuna 
        dfprod15RM[index_comuna,4] => cumsum(dfprod15RM[index_comuna,6:end]) for index_comuna in 1:53
    )
    lastepiday = semana_to_lastday.(1:77) 
    prod15map, lastepiday 
end 


interp_prod15(comuna, confmap, lastepiday) = interpolate_data(TimeArray(lastepiday, confmap[comuna])) 

grouped_prod15_confirmed(idscomunas, confmap, lastepiday) = sum(interp_prod15(id, confmap, lastepiday) for id in idscomunas)

function make_prod15_confirmed(comunas, confmap, lastepiday)
    interpolated_prod15 = Array{Float64,2}(undef, length(interp_prod15(comunas[1], confmap, lastepiday)), length(comunas))
    for i in 1:length(comunas)
        interpolated_prod15[:,i] .= interp_prod15(comunas[i], confmap, lastepiday) 
    end 
    interpolated_prod15
end 

function make_grouped_prod15_confirmed(groups, confmap, lastepiday)
    interpolated_prod15_grouped = zeros(length(grouped_prod15_confirmed(groups[1], confmap, lastepiday)), length(groups))
    for (i,group) in enumerate(groups)
        interpolated_prod15_grouped[:,i] = grouped_prod15_confirmed(group, confmap, lastepiday)
    end 
    interpolated_prod15_grouped	
end 

 
#=
function mobility_from_grouped(groupedarray::Vector{GroupData})
    number_of_groups = length(groupedarray) 
    number_of_mob_measurements = length(groupedarray[1].mobilidad) 
    mob_in_time = Array{Float64, 2}(undef, number_of_groups, number_of_mob_measurements)
    for (i, group) in enumerate(groupedarray)
        mob_in_time[i,:] .= group.mobilidad 
    end 
    mob_in_time
end 

=#


#observaciones = Array{Float64, 2}(undef, length(interpolated_data), length(comunas2) + n)
#for (i,comuna) in enumerate(comunas2)
#    observaciones[:,i] .= interpolate_data(confirmedmap[comuna])
#end 

#=
measured_index = get_index_measured_dates(timestamp(confirmedmap[comunas2[1]]))
full_grid = 1:length(measured_index)
grid = full_grid[measured_index] 

scatter!(grid, values(confirmedmap[comunas2[1]]))
=#