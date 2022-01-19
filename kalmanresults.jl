#include("KalmanMTK.jl")



results, xs, Ps = kalman_iteration(initial_u0(1.5 * a₀, gamma_e, gamma_i, beta_exterior), p);
results, xs, Ps = kalman_iteration(initial_u0(a₀, gamma_e_real, gamma_i_real, beta_exterior_real), p_real);
#results = kalman_iteration(initial_u0(a₀), p);
#@run kalman_iteration(initial_u0(a₀), p);
#include("plotting.jl")

synth_folder = is_synthetic ? "synth/" : ""
format = ".pdf"

plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown))
savefig(folder * synth_folder * "kalman_grouped_allstates_allgroups" * make_img_name(p) * format)
for class in 1:n 
    plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = class)
    savefig(folder * synth_folder * "kalman_grouped_allstates_group$class" * make_img_name(p) * ".pdf")
end 

#highlight = 1
#a_plot = plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = highlight);
#savefig(folder * synth_folder * "kalman_grouped_state$highlight" * make_img_name(p) * format)
#display(a_plot)

#for compartment in 1:5
#    plot!(a_plot, ts, [solt[n * (compartment-1) + highlight] for solt in sol.(ts)]/10^(scaling[compartment]), subplot = compartment, label = "Real solution")
#end
#plot!(ts, controls.(ts, highlight)/10^scaling[6], subplot = 6, label = "real control");
#display(a_plot)

#plot(controls.(ts, 1))
#plot!(controls.(ts, 2))

#plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = 3)
#plot_all_states_grid(tsdate, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = 5)
#@enter kalman_iteration(initial_u0(a₀), p);



#=
for compartment in 1:5
    a_plot = plot_compartment_and_incidence(tsdate, xs, Ps, states(simple_episys_uknown), total,compartment, highlight = false)

    scaling = [5., 4., 4., 5., 5., -2.];
    plot!(a_plot, tsdate, [solt[n * (compartment-1) + 1] for solt in sol.(ts)]/10^(scaling[compartment]), subplot = 1, label = "Real solution 1", legend = :topright)
    plot!(a_plot, tsdate, [solt[n * (compartment-1) + 2] for solt in sol.(ts)]/10^(scaling[compartment]), subplot = 1, label = "Real solution 2")
    plot!(a_plot, tsdate, [solt[n * (compartment-1) + 1]/total[1] for solt in sol.(ts)], subplot = 2, label = "Real solution 1")
    plot!(a_plot, tsdate, [solt[n * (compartment-1) + 2]/total[2] for solt in sol.(ts)], subplot = 2, label = "Real solution 2")
    savefig(folder * synth_folder * "kalman_grouped_$(compartment)" * make_img_name(p) * format)
end
=# 
#plot(tsdate, sol)

for state in 1:5 
    for class_to_highlight in 1:n
        a_plot = plot_compartment_and_incidence(tsdate, xs, Ps, states(simple_episys_uknown), total, state, highlight = true, class_to_highlight = class_to_highlight, estimado = is_synthetic)
        if is_synthetic
            plot_compartment_and_incidence_solution!(a_plot, tsdate, sol,
                total, state, class_to_highlight
            )
        end
        savefig(folder * synth_folder * "kalman_grouped_$(statesmap[state])_high$class_to_highlight" * make_img_name(p) * format)
    end
end


for state in 1:5 
    plot_compartment_and_incidence(tsdate, xs, Ps, states(simple_episys_uknown), total, state, highlight = false)
    savefig(folder * synth_folder * "kalman_grouped_$(statesmap[state])_allclass" * make_img_name(p) * format)
end


plot_alpha(tsdate, xs, Ps, states(simple_episys_uknown), total, highlight = false)
savefig(folder * synth_folder * "kalman_grouped_alpha_allclass" * make_img_name(p) * format)
for class_to_highlight in 1:n
    a_plot = plot_alpha(tsdate, xs, Ps, states(simple_episys_uknown), total, highlight = true, class_to_highlight = class_to_highlight, estimado = is_synthetic)
    if is_synthetic
        plot_real_control!(a_plot, tsdate, controls, class_to_highlight)
    end
    savefig(folder * synth_folder * "kalman_grouped_alpha_high$class_to_highlight" * make_img_name(p) * format)
end
#plot_alpha(tsdate, xs, Ps, states(simple_episys_uknown), total, highlight = true, class_to_highlight = 1)
#controls(0.,2)

begin 
    # gamma e 
    plot_inverted(tsdate, xs, Ps, 6n+1,
        ylabel = L"\textrm{d\acute{\imath}as}",
        label = is_synthetic ? L"\textrm{estimado}" : :none, 
        color = n+1,
        fillcolor = n+1
    ) #, title = L"1\backslash \gamma_E")
    #plot(ts, 1 ./xs[:,6n+1], title = "", ylabel = "días", label = :none)
    if is_synthetic
        hline!([1/gamma_e_real], label = L"\textrm{real}", color = n+2)
    end
    savefig(folder * synth_folder * "kalman_grouped_gamma_E_" * make_img_name(p) * format)

    # gamma i 
    plot_inverted(tsdate, xs, Ps, 6n+2,
        ylabel = L"\textrm{d\acute{\imath}as}",
        label = is_synthetic ? L"\textrm{estimado}" : :none, 
        color = n+3,
        fillcolor = n+3
    ) #, title = L"1\backslash \gamma_I")
    #savefig(folder * synth_folder * "kalman_grouped_gamma_I" * make_img_name(p) * format)
    #hline!([1/gamma_i_real], label = "γᵢ real", color = n+4)
    if is_synthetic
        hline!([1/gamma_i_real], label = L"\textrm{real}", color = n+4)
    end
    savefig(folder * synth_folder * "kalman_grouped_gamma_I_" * make_img_name(p) * format)

    # beta ext 
    plot_beta(tsdate, xs, Ps,
        label = is_synthetic ? L"\textrm{estimado}" : :none,
        color = n+5, 
        fillcolor = n+5)
    if is_synthetic
        hline!([beta_exterior_real], label = L"\textrm{real}", color = n+6) 
    end
    savefig(folder * synth_folder * "kalman_grouped_betaext_" * make_img_name(p) * format)
end
#=
b_plot = plot(title = L"\gamma_e");
plot_smoothed!(b_plot, ts, xs, Ps, states(simple_episys_uknown),
                6n+1,
                scaling_factor = 1,
                color = n+1, 
                fillcolor = n+1, 
                label = "γₑ estimado por Kalman"
            );
hline!([gamma_e_real], label = "γₑ real", color = n+4)
savefig(folder * synth_folder * "kalman_grouped_gamma_e" * make_img_name(p) * format)
display(b_plot)

b_plot = plot(title = L"\gamma_i");
plot_smoothed!(b_plot, ts, xs, Ps, states(simple_episys_uknown),
                6n+2,
                scaling_factor = 1,
                color = n+3, 
                fillcolor = n+3, 
                label = "γᵢ estimado por Kalman"
            );
hline!([gamma_i_real], label = "γᵢ real", color = n+4)
display(b_plot)


b_plot = plot(title = L"\beta_\textrm{exterior}");
plot_smoothed!(b_plot, ts, xs, Ps, states(simple_episys_uknown),
                6n+3,
                scaling_factor = 1,
                color = n+5, 
                fillcolor = n+5, 
                label = :none #"β₂ estimado por Kalman"
            );
latexify_ticks!(b_plot, ts[1], xs[1,1], [1])
hline!([beta_exterior_real], label = L"\textrm{real}", color = n+6) 
display(b_plot)
savefig(folder * synth_folder * "kalman_grouped_betaext" * make_img_name(p) * format)
=#


using Plots: hline!

#dt(odeupdater::KalmanFilter.ODEForecaster) = odeupdater.dt

plot(ts, xs[:,6n+1])

plot(results, ts, 3)

# revisar que se conserva constante la población: 
# funciona, no cambia ni en 1 persona.
plot(results.analysis[:,1:n] + results.analysis[:,n+1:2n] + results.analysis[:,2n+1:3n] + results.analysis[:,3n+1:4n])
plot(results.analysis[:,1] + results.analysis[:,n+1] + results.analysis[:,2n+1] + results.analysis[:,3n+1])
plot(results.analysis[:,n] + results.analysis[:,2n] + results.analysis[:,3n] + results.analysis[:,4n])

plot(results.analysis[:,1:n], )

plot(results.analysis[:,n+1:2n])
plot!(xs[:,n+1:2n])

plot(results.analysis[:,2n+1:3n])

plot(results.analysis[:,3n+1:4n])

plot(results.analysis[:,4n+1:5n])
plot!(xs[:,4n+1:5n])

plot(results.analysis[:,5n+1:6n], label = ["kalman 1" "kalman 2"])
plot(xs[:,5n+1:6n], label = ["smooth 1" "smooth 2"])
plot!(controls.(ts, 1), label = "real 1")
plot!(controls.(ts, 2), label = "real 2")

plot(results.analysis[:,6n+1], title = "γₑ")
plot!(xs[:,6n+1])
hline!([gamma_e_real], label = "real")

plot(results.analysis[:,6n+2], title = "γᵢ")
plot!(xs[:,6n+2])
hline!([gamma_i_real], label = "real")

plot(results.analysis[:,6n+3], title = "β₂")
plot!(xs[:,6n+3])
hline!([beta_exterior_real], label = "real")
#plot(interpolated_prod15_grouped .+ 1 , yscale = :log10)
#plot(interpolated_prod15 .+ 1, alpha = 0.2, yscale = :log10)
#plot!(observaciones .+ 1, color = [1 2 3 4 5])



plot(results.analysis[:,1])

rango = 1:length(ts)

ipss = [infocomunas[comuna].ips for comuna in comunas2] 
nombres = [infocomunas[comuna].nombre for comuna in comunas2 ]

#@run plot_all_states_grid(ts, xs, Ps, states(simple_episys_uknown), highlight = true, class_to_highlight = 3)

#savefig(folder * synth_folder * "kalman_allmuni_allstates" * make_img_name(p) * format)

b_plot = plot(Date(2020,3,4):Day(1):Date(2020,3,6), [2.,3.,1.])
latexify_ticks!(b_plot[1], :x, Date(2020,3,4))
display(a_plot)

p = plot(0:2,0:2,
    xlab = L"x^2",
    xticks = ([0,1,2],[L"3-\textrm{mar}", L"1", L"2"])
)



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

inv([
    0 1 1 1;
    1 0 1 1; 
    1 1 0 1;
    1 1 1 0 
]) * [14, 16, 24, 12]

