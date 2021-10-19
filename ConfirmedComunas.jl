#=
Datos comuna 
C:\\Users\\Tabita\\Documents\\Covid\\Datos-COVID19-MINSAL\\output\\producto15\\FechaInicioSintomas.csv
=#

using CSV 
using TimeSeries
using DataFrames 


if Sys.islinux()
    ipsfile = "data/IPS2019_por_comuna_RM.csv"
    confirmedfile = "../Datos-COVID19-MINSAL/output/producto6/bulk/data.csv" 
else
    ipsfile = "data\\IPS2019_por_comuna_RM.csv" 
    confirmedfile = "..\\Datos-COVID19-MINSAL\\output\\producto6\\bulk\\data.csv" 
end 


dfconfirmed = DataFrame(CSV.File(confirmedfile))

names(dfconfirmed)


dfconfirmedRM = filter(row -> ismissing(row["Region ID"]) ? false : row["Region ID"] .== 13. , dfconfirmed)

comunasRM = Int.(unique(dfconfirmedRM[!,"Comuna ID"]))

dfIPS = DataFrame(CSV.File(ipsfile, header = ["id", "ips"]))

struct InfoComuna
    nombre::String 
    poblacion::Float64
    codigo::Int 
    ips::Float64
end 

function getComunaInfo(comunaID::Int)
    firstdatarow = first(filter(row -> row["Comuna ID"] == comunaID, dfconfirmedRM))
    poblacion = firstdatarow[1]
    comunaID = firstdatarow[8]
    nombrecomuna = firstdatarow[9] 
    ips = filter(row -> row["id"] == comunaID, dfIPS).ips[1]
    InfoComuna(nombrecomuna, poblacion, comunaID, ips)
end 

infocomunas = Dict(
    comunaID => getComunaInfo(comunaID) for comunaID in comunasRM
) 

# comunaID => Timeseries  

splitdate(strdate) = tuple(parse.(Int, split(strdate, "/"))...)

function timeparser(strdate) 
    y, m, d = splitdate(strdate)
    Date(y, m, d)
end 

#=
function get_cases_TS(comunaID)
    casesperdate = sort(filter(row -> row["Comuna ID"] == comunaID, dfconfirmedRM), [:Fecha])[:, ["Fecha", "Casos Confirmados"]]
    TimeArray(timeparser.(casesperdate[!,1]), parse.(Float64,casesperdate[!,2]))
end 


get_cases_TS(comunasRM[1])

confirmedmap = Dict(
    comunaID => get_cases_TS(comunaID) for comunaID in comunasRM  
)  
=#

# todas estas sumas deberían dar la misma cantidad de kays de datamap (51)
#sum(length(confirmedmap[key]) .== 140 for key in keys(datamap)) # 51 
#sum(timestamp(confirmedmap[key][1])[1] .== Date(2020,03,30) for key in keys(datamap)) # 51 
#sum(timestamp(confirmedmap[key][end])[1] .== Date(2021,07,23) for key in keys(datamap)) # 51 

#==
using Plots: plot, plot! 
a_plot = plot(title = "Confirmados/población por comuna", legend = :topleft) 
[plot!(a_plot, confirmedmap[comuna]./infocomunas[comuna].poblacion, label = infocomunas[comuna].nombre) for comuna in comunas[15:20]];
a_plot  
==#