#=
Interpolate confirmed cases to all days 

Requiere haber corrido ConfirmedComunas.jl,
and confirmedmap be defined. 
=#

using Interpolations: interpolate, Gridded, Linear 
using TimeSeries

#=
Ejemplo de interpolación sencilla, que funciona con grilla
no equiespaciada. 
A = rand(20)
A_x = 1.0:2.0:40.0
nodes = (A_x,)
itp = interpolate(nodes, A, Gridded(Linear()))
scatter(A_x, A)
plot!(itp.(1.:1.:39.))
=# 

"""
    get_index_measured_dates
# Argumentos  
- `measured::Vector{Date}`: a series of no 
"""
function get_index_measured_dates(measured::Vector{Date})
    dates = first(measured):Day(1):last(measured) 
    index_measured = BitArray(zeros(length(dates)))
    for (i, date) in enumerate(dates)
        if date in measured 
            index_measured[i] = 1 
        end 
    end
    index_measured
end

#  datacomuna = confirmedmap[comunas[1]]

"""
    interpolate_data
Recibe una serie de tiempo medida cada ciertos días. Hace una estimación 
lineal de los datos en los días faltantes.
`data_comuna::TimeArray`
"""
function interpolate_data(data_comuna::TimeArray)
    measured_index = get_index_measured_dates(timestamp(data_comuna))
    full_grid = 1:length(measured_index)
    grid = full_grid[measured_index]
    data = values(data_comuna)
    itp = interpolate((grid,), data, Gridded(Linear()))
    itp.(full_grid)
end 

plot(interpolate_data(confirmedmap[comunas[1]]))
interpolate_data(confirmedmap[comunas[2]])
# quiero todos los índices de las mediciones 

A = rand(10)

scatter(measured, A) 
interpolate((measured,), A, Gridded(Linear()))

get_index_measured_dates() 