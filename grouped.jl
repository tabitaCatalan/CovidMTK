allinfos = [last(entry) for entry in infocomunas]
# encontrar y  borrar la comuna de san pedro, no hay info de mobilidad 
#for i in 1:length(allinfos) if allinfos[i].codigo == 13505 print(i) end end # 23 
allinfos = allinfos[[1:22;24:52]]

sortedbyipscomunas = sort(allinfos, by = comuna -> comuna.ips) 


ips_lims = [6, 18, 32, 45] 
ranges = ranges_from_tuple(ips_lims) 

#using Plots: hline!, vline!
#scatter(get_ips.(sortedbyipscomunas), legend = :bottomright, ms = 2, msw = 0)
#hline!([37.36, 64.37, 71.36, 77.39])
#vline!(ips_lims)
get_group_id(i) = get_id.(sortedbyipscomunas[ranges[i]])
groups = get_group_id.(1:5)

alto = GroupData("Alto", "alto", groups[1], datamap)
medioalto = GroupData("Medio Alto", "medioalto", groups[2], datamap)
medio = GroupData("Medio", "medio", groups[3], datamap)
mediobajo = GroupData("Medio Bajo", "mediobajo", groups[4], datamap)
bajo = GroupData("Bajo", "bajo", groups[5], datamap)
datagroups = [alto, medioalto, medio, mediobajo, bajo]
totales_por_clase = poblacion.(datagroups)
Pt = makePmatrix(initial_mob, mobility_from_grouped(datagroups));

# esto funciona como debería 
#plot(alto.mobilidad)
#plot!(medioalto.mobilidad)
#plot!(bajo.mobilidad)

plot(mobility_from_grouped(datagroups)')
