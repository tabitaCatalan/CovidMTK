#=
Script para la obtención de resultados con distintos casos hipotéticos 
Se necesita haber corrido el filtro de kalman con smoother, de forma que la variable xs esté definida
=#

plot(xs[:,26:30])

index_from_time(t) = floor(Int, t) + 1
# asserts 
@assert index_from_time(499.) == 500 ""
@assert index_from_time(0.) == 1 ""
@assert index_from_time(1.)== 2 ""


# comportamientos (uso de mascarilla, mantener distancia, etc)
comportamiento_normal(t, i) = xs[index_from_time(t), 5n + i]
comportamiento_cuidadoso(t, i) = xs[index_from_time(t), 5n + 1] # todos como la clase alta 
comportamiento_descuidado(t, i) = xs[index_from_time(t), 5n + 5] # todos como la clase baja

@variables i 
@register comportamiento_normal(t, i)
@register comportamiento_cuidadoso(t, i)
@register comportamiento_descuidado(t, i)

# parametros 
gamma_e_(t) = xs[index_from_time(t), 6n + 1]
gamma_i_(t) = xs[index_from_time(t), 6n + 2]
beta_2_(t) = xs[index_from_time(t), 6n + 3]
@register gamma_e_(t)
@register gamma_i_(t)
@register beta_2_(t)

# cuarentenas 
time_adjustement = 28.

function cuarentena_fuerte(t, i, j) # siguiendo al grupo que más redujo su movilidad...
    hometime = maximum([residence_times_matrix(t + time_adjustement, i, 1) for i in 1:n]) # maximo posible para ese tiempo
    if j == 1
        hometime 
    elseif j == 2
        1 - hometime
    end
end 
#cuarentena_fuerte(t,i,j) = residence_times_matrix(t + time_adjustement, 1, j) #residence_times_matrix(80., 1, j) # la cuarentena de la clase alta 
cuarentena_normal(t,i,j) = residence_times_matrix(t+time_adjustement, i, j) 
cuarentena_inexistente(t,i,j) = residence_times_matrix(0., i, j)
@register cuarentena_fuerte(t,i,j)
@register cuarentena_normal(t,i,j)
@register cuarentena_inexistente(t,i,j)



plot(ts, cuarentena_inexistente.(ts, 2, 1), label = "inexistente", title = "tiempo en casa")
for i in 1:5
    plot!(ts, cuarentena_normal.(ts, i, 1), label = "normal grupo $i")
end 
plot!(ts, cuarentena_fuerte.(ts, 2, 1), label = "fuerte")

nombres_cuarentenas = Dict(0 => "cuarentena_normal", 1 => "cuarentena_fuerte", -1 => "cuarentena_nula")
nombres_comportamientos = Dict(0 => "cuidado_normal", 1 => "cuidado_extra", -1 => "descuido")
funciones_cuarentena = Dict(0 => cuarentena_normal, 1 => cuarentena_fuerte, -1 => cuarentena_inexistente)
funciones_comportamiento = Dict(0 => comportamiento_normal, 1 => comportamiento_cuidadoso, -1 => comportamiento_descuidado)

### Simplificación estructural + resolver la EDO 

"""
# Argumentos 
- `cuarentena`:
    - `0`: normal, la gente se mueve como lo había estado haciendo
    - `-1`: cuarentena nula, no hay reducción de movilidad 
    - `1`: cuarentena fuerte, todos se mueven como el grupo 1 (movilidad mucho más reducida en general)
- `comportamiento`:
    - `0`: cuidado normal, la gente se comporta como siempre
    - `-1`: descuidado, la gente se comporta como la clase baja (que fue la que peor índice tuvo)
    - `-1`: cuidado extra, la gente se comporta como la clase alta (que fue la que mejor índice tuvo)
# Resultado 
Devuelve un diccionario con keys `"nombre"`, `"cuarentena"`, `"comportamiento"`.
"""
function make_cases_metadata(tipo_cuarentena, tipo_comportamiento)
    
    name = string(tipo_cuarentena) * "_"  * string(tipo_comportamiento) * "_" * nombres_cuarentenas[tipo_cuarentena] * "_y_" * nombres_comportamientos[tipo_comportamiento]
    Dict(
        "nombre" => name,
        "cuarentena" => funciones_cuarentena[tipo_cuarentena],
        "comportamiento" => funciones_comportamiento[tipo_comportamiento]
    )
end 

#casos = [(0,0), (0,1), (0,-1), (-1, 0), (-1, 1), (-1, -1), (1, 0), (1, 1), (1, -1)]
casos = [
    (-1, -1), (-1, 0), (-1, 1), 
    (0,-1), (0,0), (0,1), 
    (1, -1), (1, 0), (1, 1)
]

metadata_casos = [make_cases_metadata(cuarentena, comportamiento) for (cuarentena, comportamiento) in casos]

function make_and_solve_case(t, n, m, metadata_caso, u0, p)
    @named modelo = epi_model_known_all_inputs(t, n, m, metadata_caso["cuarentena"], metadata_caso["comportamiento"], gamma_e_, gamma_i_, beta_2_)
    modelo_simple = structural_simplify(modelo);
    problema = ODEProblem(modelo_simple, u0, (0., T-dt), p)
    sol = solve(problema, Tsit5(), saveat = 1.);
    Dict(
        "modelo" => modelo,
        "modelo_simple" => modelo_simple,
        "problema" => problema, 
        "solucion" => sol
    )
end

# posiblemente cambiar C0 por 0s 
final_u0 = make_x_k(episys_uknown, xs[1,1:5], xs[1,6:10], xs[1,11:15], xs[1,16:20], xs[1,21:25])
final_p = p

solutions = [make_and_solve_case(t, n, m, metadata, final_u0, final_p) for metadata in metadata_casos];


a_plot = plot(framestyle=:box, size = (1800, 1200), layout = (3,3), link = :y);
for i in 1:length(casos)
    plot!(a_plot, solutions[i]["solucion"],
        vars = (t, E),
        title = "Infectados " * metadata_casos[i]["nombre"], 
        subplot = i, label = :none, ylims = (0., 3e4), ylabel = "")
end
latexify_ticks!(a_plot, ts[1], " ")
plot!(a_plot, title = ["Descuido" "Cuidado normal" "Cuidado extra" "" "" "" "" "" ""], titlefont = font(20))
plot!(a_plot, ylabel = ["Sin cuarentena" "" "" "Cuarentena normal" "" "" "Cuarentena extra" "" ""], left_margin = 10mm, yguidefontsize=20)
display(a_plot)
savefig(folder * "hipcases" * make_img_name(p) * ".pdf")

#latexify_ticks!(a_plot[1], :x, 1.)
# De estos mencionar solamente la proporción final, el máximo y el tiempo estimado que tarda 
plot(solutions[1]["solucion"], vars = (t, S)) # la población de enferma dentro de los primeros 3 meses hasta alcanzar la inmunidad de rebaño 
plot(solutions[2]["solucion"], vars = (t, S)) 


incidence(t, E) = (t, E./total)

toromanlatex(str) =to_latex_string("\\mathrm{" * replace(str, " " => "\\,\\,") * "}")

function plot_all_hipot_cases(solutions, ts, index; xymeasures = (450, 200), titlefont = 12, normalize = false, common_ylims = false, ylims = (0.,1.))
    xmeasure, ymeasure = xymeasures
    
    #a_plot = plot(framestyle=:box, size = (1800, 1200), layout = (3,3), link = :y);
    a_plot = plot(framestyle=:box, size = (xmeasure * 3, ymeasure * 3), layout = (3,3), link = :y);
    for i in 1:length(casos)
        data = (solutions[i]["solucion"][index, :])';
        if normalize
            data ./= total' 
        end
        plot!(a_plot,
            ts,
            data,
            #title = "Incidencia Infectados " * metadata_casos[i]["nombre"], 
            subplot = i, label = :none, ylabel = "")
        if common_ylims
            plot!(a_plot, ylims = ylims)
        end
    end
    latexify_ticks!(a_plot, ts[1], " ")
    plot!(a_plot, title = toromanlatex.(["Descuido" "Cuidado normal" "Cuidado extra" "" "" "" "" "" ""]), titlefont = font(titlefont))
    plot!(a_plot, ylabel = toromanlatex.(["Sin cuarentena" "" "" "Cuarentena normal" "" "" "Cuarentena extra" "" ""]), left_margin = 10mm, yguidefontsize=titlefont)
    a_plot
end



function plot_compare_with_normal_case(solutions, ts, index, case; xymeasures = (450, 200), normalize = false)
    normal_case = 5
    #a_plot = plot(framestyle=:box, size = (1800, 1200), layout = (3,3), link = :y);
    a_plot = plot(framestyle=:box, size = xymeasures);
    data_normal = (solutions[normal_case]["solucion"][index, :])';
    if normalize
        data ./= total' 
    end    
    data_case = (solutions[case]["solucion"][index, :])';
    if normalize
        data_case ./= total' 
    end    
    
    plot!(a_plot,
        ts,
        data_normal,
        label = :none, ylabel = "", alpha = 0.3)
    plot!(a_plot,
        ts,
        data_case,
        label = :none, ylabel = "", color = [1 2 3 4 5])
    latexify_ticks!(a_plot, ts[1], " ")
    
    a_plot
end


indexI = 11:15
a_plot = plot_all_hipot_cases(solutions, tsdate, indexS; normalize = true)
plot_all_hipot_cases(solutions, tsdate, indexS; normalize = true)
savefig(folder * "allhipcases_S-N_" * make_img_name(p) * ".pdf")

plot_all_hipot_cases(solutions, tsdate, indexI; normalize = true, common_ylims = true, ylims = (0., 0.15))
savefig(folder * "allhipcases_I-N_commonylim0-15_" * make_img_name(p) * ".pdf")

plot_all_hipot_cases(solutions, tsdate, indexI; normalize = true, common_ylims = false)
savefig(folder * "allhipcases_I-N_" * make_img_name(p) * ".pdf")

plot_all_hipot_cases(solutions, tsdate, indexI; normalize = true, common_ylims = true, ylims = (0.,0.015))
savefig(folder * "allhipcases_I-N_commonylim0-015" * make_img_name(p) * ".pdf")

plot_all_hipot_cases(solutions, tsdate, indexI; normalize = false, common_ylims = true, ylims = (0., 5e4))
savefig(folder * "allhipcases_I_commonylim5e4_" * make_img_name(p) * ".pdf")



function plot_important_dates!(a_plot, important_dates; relx = 0.03, rely = 0.9, fontsize = 8, letras = true)
    char = letras ? 'a' : 1
    for i in 1:length(important_dates)
        vline!(a_plot, important_dates[i:i], label = :none, c = :grey)
        relxbase = to_relative_x(a_plot[1], important_dates[i])
        if letras 
            texto = to_latex_string("\\textrm{($(char))}")
        else 
            texto = to_latex_string("($(char))")
        end 
        annotate!(a_plot, (relxbase + relx, rely), text(texto, fontsize))
        char += 1
    end
    a_plot
end

# necesito la transformación inversa de esta 
# a partir del valor absoluto, necesito el relativo 
# Ideas de 
# https://discourse.julialang.org/t/plots-xlims-with-x-axis-given-by-dates/47926
# https://github.com/JuliaPlots/Plots.jl/issues/2728
function to_relative_x(p::Plots.Subplot, datex)
    x = Dates.value(datex)
    xlims = Plots.xlims(p)
    (x - xlims[1])/(xlims[2] - xlims[1])
end 

important_dates = [Date(2020, 5, 5), Date(2021, 3, 1), Date(2021, 6, 9)]
for i in [3, 4, 7, 8] # Casos interesantes
    a_plot = plot_compare_with_normal_case(solutions, tsdate, indexI, i)
    plot_important_dates!(a_plot, important_dates; letras = false)

    savefig(folder * "comparecase_$(i)withnormal_I_" * make_img_name(p) * ".pdf")
end
