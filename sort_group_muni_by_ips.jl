
calculate_population(ids) = sum(infocomunas[id].poblacion for id in ids)


#=
function find_opt_limranges()
    minsol = Inf; soltuple = (1,1)
    for i in 10:40 
        for j in i+1:41
            val = sum(abs.(population_by_class(ranges_from_tuple((i,j)), sortedbyipscomunas) .-  totalpopulation/social_classes)) 
            if val < minsol
                minsol = val 
                soltuple = (i,j)
                println("Val = ", val)
            end 
        end
    end 
end 
=# 


"""
    ranges_from_tuple(lims) 
devuelve los rangos de los distintos grupos a partir de los índices de corte.
Supone que lims es una tupla de elementos ordenados en forma ascendente.
- `lims`: tuple de índices de corte para la selección de comunas. 
- `last`: last index (default 51)
# Ejemplo 
```julia
julia> ranges_from_tuple((15, 30))
3-element Vector{UnitRange{Int64}}:
 1:15
 16:30
 31:51
```
"""
function ranges_from_tuple(lims, last = 51) 
    if !issorted(lims)
        println("Not sorted lims: $lims")
    end
    ranges = [1:lims[1]]
    for i in 1:length(lims)-1
        push!(ranges, lims[i]+1:lims[i+1])
    end  
    push!(ranges, lims[end]:last)
    ranges 
end 

# Necesito hacer la matriz y las condiciones iniciales
# para datos agrupados... 
# Además de hacer gráficos >.< no hay tiempo  

function grouped_mobility(idscomunas, datamap)
    total_clase = calculate_population(idscomunas)
    sum(infocomunas[id].poblacion .* datamap[id] for id in idscomunas) ./ total_clase 
end 

get_id(info) = info.codigo
get_ips(info) = info.ips



struct GroupData
    """Nombre del grupo, para usar en títulos de gráficos, etc."""
    nombre::String
    """Abreviación del nombre, para usar en labels de gráficos, etc."""
    abr::String
    """Lista de ids (códigos) de todas las municipalidades del grupo."""
    ids::Vector{Int}
    """Movilidad agrupada, calculada ponderando la movilidad de cada municipalidad por su fracción de la población dentro del grupo."""
    mobilidad::Vector{Float64}
    """Total de población del grupo"""
    poblacion::Float64
end 


function GroupData(nombre, abr, ids, datamap)
    GroupData(nombre, abr, ids, grouped_mobility(ids, datamap), calculate_population(ids))
end 

poblacion(gd::GroupData) = gd.poblacion

function mobility_from_grouped(groupedarray::Vector{GroupData})
    number_of_groups = length(groupedarray) 
    number_of_mob_measurements = length(groupedarray[1].mobilidad) 
    mob_in_time = Array{Float64, 2}(undef, number_of_groups, number_of_mob_measurements)
    for (i, group) in enumerate(groupedarray)
        mob_in_time[i,:] .= group.mobilidad 
    end 
    mob_in_time
end 
