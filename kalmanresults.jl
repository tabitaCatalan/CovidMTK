results, xs, Ps = kalman_iteration(initial_u0(a₀), p);

plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown))
plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = 4)
#@enter kalman_iteration(initial_u0(a₀), p);

plot(interpolated_prod15_grouped .+ 1 , yscale = :log10)
plot(interpolated_prod15 .+ 1, alpha = 0.2, yscale = :log10)
plot!(observaciones .+ 1, color = [1 2 3 4 5])

plot(results.analysis[:,1])

rango = 1:length(ts)


using Plots: plot, plot!, cgrad 

ipss = [infocomunas[comuna].ips for comuna in comunas2] 
nombres = [infocomunas[comuna].nombre for comuna in comunas2 ]




a_plot = plot(layout=(3,2),framestyle=:box, link = :x, size = (800, 450));
for state in 1:5 # estados 
    for i = 1:n # comunas 
        plot_smoothed!(a_plot, tsdate, xs, Ps, 1., (state-1)*n + i, subplot = state, fillalpha = 0.1) # 10^5 
    end 
end 
plot_smoothed!(a_plot, tsdate, xs, Ps, 1., 5*n + 1, subplot = 6, fillalpha = 0.1); # 10^5 
if ! one_control 
    for i = 2:n # comunas 
        plot_smoothed!(a_plot, tsdate, xs, Ps, 1., 5*n + i, subplot = 6, fillalpha = 0.1) # 10^5 
    end 
end 
display(a_plot)

#savefig(folder * "kalman_allmuni_allstates" * make_img_name(p) * ".svg")
d_plot = plot(framestyle=:box, size = (800, 450));
#for state in 1:6 # estados  

plot_smoothed!(d_plot, ts, xs, Ps, 1., 5*n + 1, fillalpha = 0.1) # 10^5 
if ! one_control 
    for i = 2:n # comunas 
        plot_smoothed!(d_plot, ts, xs, Ps, 1., 5*n + i,  fillalpha = 0.1) # 10^5 
    end 
end 
plot!(ts, results.analysis[:,26])
display(d_plot)

function plot_highlighted_muni(ts, xs, Ps, total, i0, state)
    cg = cgrad();
    extra = 0.5
    colors = cg[normalize_extra(ipss, extra)];

    e_plot = plot(framestyle=:box, size = (800, 450));
    for i = 1:n # comunas 
        plot_smoothed!(e_plot, ts, xs, Ps, total[i], (state-1)*n + i, fillalpha = 0.15, fillcolor = :gray90, color = :grey90, label = :none) # 10^5 
    end 
    plot_smoothed!(e_plot, ts, xs, Ps, total[i0], (state-1)*n + i0,
            label = infocomunas[comunas2[i0]].nombre,
            line_z = infocomunas[comunas2[i0]].ips, 
            clims = extralims(ipss, extra),
            color = colors[i0],
            fillcolor = colors[i0],
            fillalpha = 0.15
        )
    e_plot
end 

display(plot_highlighted_muni(ts, xs, Ps, total, 5, 2))


function extralims(array, out)
    minarray = minimum(array)
    maxarray = maximum(array) 
    extra = out * (maxarray - minarray) / 2 
    minarray - extra, maxarray + extra
end 

function normalize_extra(array, exlims:: NTuple{2, Float64})
    exmin, exmax = exlims 
    (array .-exmin)/(exmax - exmin)
end 

normalize_extra(array, out) = normalize_extra(array, extralims(array, out))


b_plot = plot(title = "Percentage susceptibles");
for i in 1:n
    plot_smoothed!(b_plot, ts, xs, Ps, total[i], i,
        label = infocomunas[comunas2[i]].nombre,
        line_z = infocomunas[comunas2[i]].ips, 
        clims = extralims(ipss, extra),
        color = colors[i],
        fillcolor = colors[i],
        fillalpha = 0.1
    )
end 
display(b_plot) 

plot_smoothed_with_cgrad!(b_plot, ts, xs, Ps, total, i, extra)

b_plot = plot(title = "Incidencia");
for i in 1:3
    plot_smoothed!(b_plot, ts, xs, Ps, totales_por_clase[i], n + i,
        fillalpha = 0.1
    )
end 
display(b_plot) 


