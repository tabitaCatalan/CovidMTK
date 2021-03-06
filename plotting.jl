#=
Aux functions for plotting 
=#
using Plots
using Plots.Measures: mm
using LaTeXStrings
using Printf: @sprintf 
#======================================
Working with spanish days and months 
https://docs.julialang.org/en/v1/stdlib/Dates/
See locale attribute 
======================================#

spanish_months = ["enero", "febrero", "marzo", "abril", "mayo", "junio", "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre"];
spanish_months_abbrev = ["ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic"];
spanish_days = ["lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo"];
spanish_days_abbrev = ["lun", "mar", "mie", "jue", "vie", "sáb", "dom"];

Dates.LOCALES["spanish"] = Dates.DateLocale(spanish_months, spanish_months_abbrev, spanish_days, spanish_days_abbrev);

#=
Aux funtion to treat MTK variables 
=#
"""
Return a Tuple with string of a Symbolic variable name and a string of its index.
# Example 
```
julia> @variables t R[1,2](t)
julia> var = R[2]
R[2](t)
julia> split_var_and_real_index(var)
("R", "2")
```
"""
function split_var_and_real_index(var::Num)
    longstrvar = collect(string(var))
    leftbracketindex = findfirst(isequal('['), longstrvar)
    rightbracketindex = findfirst(isequal(']'), longstrvar)
    strvar = join(longstrvar[1:leftbracketindex-1])
    strindex = join(longstrvar[leftbracketindex+1:rightbracketindex-1])
    strvar, strindex
end 

#=============================== #
# -------------------------------#
# Symbolics variables to Strings #
# -------------------------------#
# ===============================# 

#=============================== #
Subindex strings 
# ===============================# 
subindexs = Dict(
    "0" => "₀",
    "1" => "₁",
    "2" => "₂",
    "3" => "₃",
    "4" => "₄",
    "5" => "₅",
    "6" => "₆",
    "7" => "₇",
    "8" => "₈",
    "9" => "₉",
)


"""
Transform a string of a number to a string of the number using subindexs.
# Example 
```julia
julia> to_subindex("1729")
"₁₇₂₉"
```
"""
to_subindex(strnum) = join([subindexs[digit] for digit in split(strnum, "")])


"""
Return a string from a symbolic ModelingToolkit array variable using subindexs.
- `var`: symbolic array ModelingToolkit variable 
- `index`: index number 
# Example 
```julia
@variables t S[1,2](t)
string(S[1]) #"S[1](t)"
to_string(S[1]) #"S₁(t)"
```
"""
function to_subindex_string(var::Num)
    strvar, strindex = split_var_and_real_index(var)
    strvar * to_subindex(strindex) * "(t)"
end  

#====================================== #
Latex strings 
# ======================================#

"""
Return a Latex string from a symbolic ModelingToolkit array variable.
- `var`: symbolic ModelingToolkit variable 
- `index`: index number 
```
"""
function to_latex_string(var::Num)
    strvar, strindex = split_var_and_real_index(var)
    strvar = replace(strvar, "α" => "\\alpha")
    L"%$(strvar)_{%$(strindex)}(t)"
end  

#=============================== #
# -------------------------------#
# Latexify                       #
# -------------------------------#
# ===============================# 

#=============================== #
 Numbers and dates 
# ===============================#

function to_latex_string(t::Date)
    (y, m, d) = Dates.yearmonthday(t)
    strm = Dates.monthabbr(m; locale = "spanish")
    L"%$d \textrm{/%$(strm)/} %$y"
end 

to_latex_string(ft::Float64) = latexstring(@sprintf "%.2f" ft)

function to_latex_string(a)
    latexstring(a)
end 

#=============================== #
 Ticks in plots 
# ===============================#

# El tipo de `example_data` define el comportamiento de `transform`
transform(strdate, example_data::Date) = Date(strdate, "y-m-d")
transform(strfloat, example_data::Float64) = parse(Float64, replace(strfloat, "−" => "-"))
transform(strfloat, example_data::String) = replace(replace(strfloat, "−" => "-"), "×" => "\\times")
"""
Replace ticks from x or y axis of a simple Plots with a LaTeXString version.
# Arguments 
- `a_plot::Plots.SubPlot`: a simple Plot, this not intended to work with plots with subplots.
- `axis::Symbol`: options are `:x` and `:y`. 
- `y0`
"""
function latexify_ticks!(a_plot::Plots.Subplot, axis::Symbol, example_data)
    if axis == :x  
        oldticks = Plots.xticks(a_plot)
    elseif axis == :y
        oldticks = Plots.yticks(a_plot)
    end 
    newticks = to_latex_string.(transform.(oldticks[2], example_data))
    if axis == :x  
        plot!(a_plot, xticks = (oldticks[1], newticks))
    elseif axis == :y
        plot!(a_plot, yticks = (oldticks[1], newticks))
    end 
end 

"""
Applies `latexify_ticks!` in x and y axis to every subplot of `a_plot`
"""
function latexify_ticks!(a_plot::Plots.Plot, t0, y0, notremoved = 1:length(a_plot))
    for subplot in 1:length(a_plot)
        latexify_ticks!(a_plot[subplot], :y, y0)
    end 
    for subplot in notremoved
        latexify_ticks!(a_plot[subplot], :x, t0)
    end
end
#===========================
Plot with scientific notation:
===========================# 
#=
# Conservé este código porque es útil para graficar el caso sintético.

calculate_exponent(ymin, ymax) = round(Int, log10(ymax - ymin)) + 1  
calculate_exponent(ymin, ymax) = floor(Int, log10(ymax)) 
minmax(sol, var) = (minimum([sol[var[1]] sol[var[2]]]), maximum([sol[var[1]] sol[var[2]]]))

function plot_scnotation!(a_plot, sol, var, subplot, scat = false)
    a = calculate_exponent(minmax(sol,var)...)
    for index in 1:length(var)
        label = to_string(var, index)
        if scat
            scatter!(a_plot, sol.t, sol[var[index]]/10^a, label = label, ylabel = "people", subplot = subplot, xticks = :none, ms = 1, msw = 0)
        else 
            plot!(a_plot, sol.t, sol[var[index]]/10^a, label = label, ylabel = "people", subplot = subplot, xticks = :none)
        end 
    end 
    #annotate!(a_plot, (-0.085, 0.98), text("x10^$a", :left, 9)) # esto funcionaba con 1 plot 
    annotate!(a_plot, (-0.08, 1.1), text("x10^$a", :left, 8), subplot = subplot)    
end 
=#

function add_scix10_to_plot!(a_plot, subplot, exponent::Int)
    if exponent != 0
        annotate!(a_plot, (-0.07, 1.1), text(L"\times 10^{%$(exponent)}", :left, 8), subplot = subplot)
    end
end

vars_in_state(state) = ((state-1)*n + 1):((state-1)*n + n)
calculate_exponent(ymax) = floor(Int, log10(ymax)) 

function calculate_scaling_exponents(xs, state)
    ymax = maximum(xs[:,vars_in_state(state)])
    calculate_exponent(ymax)
end

#============================= # 
A versatil plotting function 
# =============================#


"""
    plot_smoothed!(a_plot, ts, xs, Ps, a, index; kwargs...)
Grafica un serie de datos con barra de error. Agrega un label con el nombre del estado, suponiendo que el sistema es `simple_episys_uknown`.
# Argumentos 
- `a_plot::Plots.Plot`, se sobreescribirá para agregar el nuevo gráfico 
- `ts`: datos del eje x 
- `xs::Array{T,2}`: datos del eje y, con tantas filas con valores en `ts`. Cada columna corresponde a un estado distinto. 
- `Ps::Array{T,3}`: matriz de varianzas covarianza del vector de estados en el tiempo.
    `Ps[:,:,n]` corresponde a la matriz de covarianzas en tiempo `ts[n]`, por lo que tiene tantas filas y columnas 
    como columnas tiene `xs`. Se usan solo las diagonales de la forma `Ps[index, index, n]` (para el error del estado `index` a tiempo `n`).
- `symstates`: array of Modeling Toolkit Symbolic Variables.
- `index`: se graficará solo el estado `index`-ésimo.
## Opcionales
- `scaling_factor = 1.`: un factor de escalamiento (útil para evitar que aparezca números con notación científica como `x 10^4`).
- `kwargs`: keyword arguments que se pasan a la función `plot!`.
"""
function plot_smoothed!(a_plot, ts, xs, Ps, symstates, index; scaling_factor = 1., kwargs...)
    #label = to_latex_string(symstates[index]) # antes era "s$index"
    plot!(ts, xs[:,index] * scaling_factor, ribbon = sqrt.(Ps[index,index,:]) * scaling_factor; kwargs...)
end 

#================================= #
Aux functions useful in plotting 
# =================================#

remove_xticks(state) = ! (state in 5:6)

remove_xticks!(a_plot, subplot) = plot!(a_plot[subplot], xticks = (Plots.xticks(a_plot[subplot])[1], ["" ]))

"""
Devuelve un Diccionario con attributos para graficar cierta clase 
# Argumentos 
- `class::Int`: clase a graficar
- `var::Num`: MTK variable, corresponiente a la clase y estado a graficar.
- `highlight::Bool`: si es `true`, entonces se supondrá que se hará un gráfico que destaca la clase `class_to_highlight`.
    Todas las demás se grafican en gris. Si es `false` cada clase toma un color asociado a su clase (tanto para los datos 
    como para la barra de error).
- `class_to_highlight`: solo se usa si `highlight` es `true`.
"""
function calculate_plot_attrib(class::Int, var::Num, highlight::Bool, class_to_highlight::Int)
    if highlight && class_to_highlight != class
        color = :grey90
        fillcolor = :gray90
        label = :none
        fillalpha = 0.3
    else 
        fillalpha = 0.15
        color = class
        fillcolor = class
        label = to_latex_string(var) # antes era "s$index"
    end
    Dict(:color => color, :fillcolor => fillcolor, :fillalpha => fillalpha, :label => label)
end



"""
    put_at_the_end(array, index)
Returns an array with the `index`-th element removed and inderted at the end.
# Example 
```
julia> put_at_the_end(1:5, 3)
5-element Vector{Int64}:
 1
 2
 4
 5
 3
```
"""
put_at_the_end(array, index) = array[[1:(index-1); (index+1):end; index]]

#================================= # 
High level plotting functions 
# =================================#

"""
    plot_all_states_grid(ts, xs, Ps)
Dibuja una grilla de (3,2) para cada uno de los compartimientos de un sistema SEIR. 
- (1,1): Susceptibles (S)
- (1,2): Expuestos (E)
- (2,1): Infectados (I)
- (2,2): Recuperados (R)
- (3,1): Acumulados (C)
- (3,2): Tasa de contagio (α)
"""
function plot_all_states_grid(ts, xs, Ps, symstates; highlight = false, class_to_highlight = n)
    class_to_highlight = highlight ? class_to_highlight : n
    scaling_exponents = [calculate_scaling_exponents(xs, state) for state in 1:6]
    scaling_factors = 10 .^(-Float64.(scaling_exponents))
    a_plot = plot(layout=(3,2),framestyle=:box, link = :x, size = (800, 450));
    if one_control
        total_states_to_plot = 5
        plot_smoothed!(a_plot, ts, xs, Ps, symstates, 5*n + 1,
            scaling_factor = scaling_factors[state],
            subplot = 6,
            fillalpha = 0.1,
            color = :black
        ); 
    else 
        total_states_to_plot = 6
    end
    for state in 1:total_states_to_plot # estados 
        for class = put_at_the_end(1:n, class_to_highlight) # clases  
            index = (state-1)*n + class
            attr = calculate_plot_attrib(class, Num(symstates[index]), highlight, class_to_highlight)
            plot_smoothed!(a_plot, ts, xs, Ps, symstates,
                index,
                scaling_factor = scaling_factors[state],
                subplot = state,
                fillalpha = attr[:fillalpha],
                color = attr[:color], 
                fillcolor = attr[:fillcolor], 
                label = attr[:label]
            ) # 10^5 
        end 

        add_scix10_to_plot!(a_plot, state, scaling_exponents[state])
        if remove_xticks(state)
            remove_xticks!(a_plot, state)
        end
    end 
    plot!(a_plot, top_margin = 3mm) 
    latexify_ticks!(a_plot, ts[1], xs[1,1], [5,6])
    a_plot 
end 


function plot_compartment_and_incidence(ts, xs, Ps, symstates, totals, state; highlight = false, class_to_highlight = n, estimado = false, title = false)
    class_to_highlight = highlight ? class_to_highlight : n
    scaling_exponents = [calculate_scaling_exponents(xs, state) for state in 1:6]
    scaling_factors = 10 .^(-Float64.(scaling_exponents))
    #a_plot = plot(layout=(1,2),framestyle=:box, link = :x, size = (800, 150)); <- este funciona
    #a_plot = plot(layout=(2,1),framestyle=:box, link = :x, size = (500, 300), palette = palette([:green, :yellow, :red], 5));
    #a_plot = plot(layout=(2,1),framestyle=:box, link = :x, size = (500, 300), palette = palette([:lightseagreen, :darkgoldenrod1, :firebrick], 5));
    a_plot = plot(layout=(2,1),framestyle=:box, link = :x, size = (500, 300));
    
    for class = put_at_the_end(1:n, class_to_highlight) # clases  
        index = (state-1)*n + class
        attr = calculate_plot_attrib(class, Num(symstates[index]), highlight, class_to_highlight)
        plot_smoothed!(a_plot, ts, xs, Ps, symstates,
            index,
            scaling_factor = 1., #scaling_factors[state],
            subplot = 1,
            fillalpha = attr[:fillalpha],
            color = attr[:color], 
            fillcolor = attr[:fillcolor], 
            label = estimado && class == class_to_highlight ? L"\textrm{estimado}" : :none,
            #legend = :none, 
            title = title ? L"\textrm{%$(statesfullnamemap[state])}" : "",
            ylabel = L"\,%$(statesmap[state])_i(t)"
        ) # 10^5 
        incidlabel = attr[:label] == :none ? :none : L"%$(statesmap[state])_{%$(class)}(t) /N_{%$(class)}" 
        
        plot_smoothed!(a_plot, ts, xs, Ps, states(simple_episys_uknown),
            index,
            scaling_factor = 1/totals[class],
            subplot = 2,
            color = attr[:color], 
            fillcolor = attr[:fillcolor], 
            fillalpha = attr[:fillalpha],
            label = estimado && class == class_to_highlight ? L"\textrm{estimado}" : "",
            #legend = :none, 
            #ylabel = L"\textrm{%$(statesfullnamemap[state])}\, %$(statesmap[state])_i(t) /N_i"
            ylabel = L"%$(statesmap[state])_i(t) /N_i"
        );
    end 

    #add_scix10_to_plot!(a_plot, 1, scaling_exponents[state])

    plot!(a_plot, top_margin = 3mm) 
    latexify_ticks!(a_plot[1], :y, " ")
    latexify_ticks!(a_plot[1], :x, ts[1])
    latexify_ticks!(a_plot[2], :y, xs[1,1])
    latexify_ticks!(a_plot[2], :x, ts[1])
    a_plot 
end

function plot_compartment_and_incidence_solution!(a_plot, ts, sol, total, state, class_to_hightlight)
    plot!(a_plot, ts, sol[n*(state - 1)+class_to_hightlight, 1:end-1],
        subplot = 1, label = L"\textrm{real}", 
        color = n+6 + class_to_hightlight)
    plot!(a_plot, ts, sol[n*(state - 1)+class_to_hightlight, 1:end-1] / total[class_to_hightlight], subplot = 2,
        label = L"\textrm{real}",
        color = n+6 + class_to_hightlight
    )
end 


function plot_alpha(ts, xs, Ps, symstates, totals; highlight = false, class_to_highlight = n, title = false, estimado = false)


    scaling_exponents = [calculate_scaling_exponents(xs, state) for state in 1:6]
    scaling_factors = 10 .^(-Float64.(scaling_exponents))
    #a_plot = plot(layout=(1,2),framestyle=:box, link = :x, size = (800, 150)); <- este funciona
    #a_plot = plot(layout=(2,1),framestyle=:box, link = :x, size = (500, 300), palette = palette([:green, :yellow, :red], 5));
    #a_plot = plot(layout=(2,1),framestyle=:box, link = :x, size = (500, 300), palette = palette([:lightseagreen, :darkgoldenrod1, :firebrick], 5));
    a_plot = plot(framestyle=:box, size = (500, 160));
    
    state = 6
    for class = put_at_the_end(1:n, class_to_highlight) # clases  
        index = (state-1)*n + class
        attr = calculate_plot_attrib(class, Num(symstates[index]), highlight, class_to_highlight)
        plot_smoothed!(a_plot, ts, xs, Ps, symstates,
            index,
            scaling_factor = 1., #scaling_factors[state],
            fillalpha = attr[:fillalpha],
            color = attr[:color], 
            fillcolor = attr[:fillcolor], 
            label = estimado && class == class_to_highlight ? L"\textrm{estimado}" : :none, 
            title = title ? L"\textrm{Factor}\,\,\textrm{sanitario}\,\,\alpha_i" : "",
            #ylabel = L"\,%$(statesmap[state])_i(t)"
        ) # 10^5 
    end
    #add_scix10_to_plot!(a_plot, 1, scaling_exponents[6])
    plot!(ylims = (0., 0.03))
    latexify_ticks!(a_plot, ts[1], xs[1,1], [1])
    a_plot
end

function plot_real_control!(a_plot, ts, controls, class_to_highlight)
    plot!(a_plot, ts, controls.(0.:dt:T-dt, class_to_highlight), color = n+6 + class_to_highlight, label = L"\textrm{real}")
end


error_propagation_add((x, Δx), (y, Δy)) = (x + y, √(Δx^2 + Δy^2))

relative_error((x, Δx), (y, Δy)) = √((Δx/x)^2 + (Δy/y)^2)

function error_propagation_mul((x, Δx), (y, Δy))
    q = xy 
    Δq = abs(q) * relative_error((x, Δx), (y, Δy))
    (q, Δq)
end
function error_propagation_div((x, Δx), (y, Δy))
    q = x/y 
    Δq = abs(q) * relative_error((x, Δx), (y, Δy))
    (q, Δq)
end

get_data_and_error(index, t, xs, Ps) = (xs[t, index], sqrt(Ps[index, index, t]))

function plot_inverted(ts, xs, Ps, index; kwargs...)

    a_plot = plot(framestyle=:box, size = (500, 160));
    gamma_and_error = [error_propagation_div((1.,0.), get_data_and_error(index, k, xs, Ps)) for k in 1:Nmediciones]
    plot!(a_plot, ts, [gamma_and_error[i][1] for i in 1:Nmediciones], ribbon = [gamma_and_error[i][2] for i in 1:Nmediciones]; kwargs...)

    
    latexify_ticks!(a_plot, ts[1], xs[1,1], [1])
    a_plot
end

function plot_beta(ts, xs, Ps; kwargs...)

    a_plot = plot(framestyle=:box, size = (500, 160));
    plot_smoothed!(b_plot, ts, xs, Ps, states(simple_episys_uknown),
                6n+3,
                scaling_factor = 1,
                color = n+5, 
                fillcolor = n+5;
                kwargs... #"β₂ estimado por Kalman"
            );
    latexify_ticks!(a_plot, ts[1], xs[1,1], [1])
    a_plot
end

statesmap = ["S", "E", "I", "R", "C", "\alpha"]
statesfullnamemap = ["Susceptibles", "Expuestos", "Infectados", "Recuperados", "Acumulados", ""]

#=s1 = scatter(1:10, 1:10, legend = false, color = 1)
s2 = scatter(1:10, 1:10, legend = false, color = 2)
s3 = scatter(1:10, 1:10, legend = false, color = 3)
s4 = scatter(1:10, 1:10, legend = false, color = 4)
legend = plot([0 0 0 0], showaxis = false, grid = false, label = ["s1" "s2" "s3" "s4"])
plot(s1, s2, s3, s4, legend, layout = @layout([[A B; C D] E{.1w}])) 
=#
