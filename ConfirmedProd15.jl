#=
confirmed using product 15 (por semana epidemiológica)
=#

using CSV 
using TimeSeries
using DataFrames 

prod15file = "C:\\Users\\Tabita\\Documents\\Covid\\Datos-COVID19-MINSAL\\output\\producto15\\FechaInicioSintomas.csv"
dfprod15 = DataFrame(CSV.File(prod15file))

dfprod15RM = filter(row -> ismissing(row["Codigo region"]) ? false : row["Codigo region"] .== 13. , dfprod15)

prod15map = Dict(
    dfprod15RM[index_comuna,4] => Array(dfprod15RM[index_comuna,6:end]) for index_comuna in 1:53
)

semana_to_lastday(week) = Date(2020,2,22) + (week - 1) * Day(7)
semana_to_lastday(1)
#plot(cumsum(prod15map[comunas2[1]]))
#plot!(cumsum(prod15map[comunas2[2]]))


interp_prod15(comuna) = interpolate_data(TimeArray(semana_to_lastday.(1:77), cumsum(prod15map[comuna]))) 

interpolated_prod15 = [interp_prod15(comunas2[1]) interp_prod15(comunas2[2])]

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