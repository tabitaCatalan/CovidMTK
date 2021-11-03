#=
Aux functions for plotting 
=#
using Plots
using LaTeXStrings
#===========================
Symbolics variables to Strings 
===========================# 
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
function split_var_and_real_index(var)
    longstrvar = collect(string(var))
    leftbracketindex = findfirst(isequal('['), longstrvar)
    rightbracketindex = findfirst(isequal(']'), longstrvar)
    strvar = join(longstrvar[1:leftbracketindex-1])
    strindex = join(longstrvar[leftbracketindex+1:rightbracketindex-1])
    strvar, strindex
end 

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
function to_subindex_string(var)
    strvar, strindex = split_var_and_real_index(var)
    strvar * to_subindex(strindex) * "(t)"
end  

"""
Return a Latex string from a symbolic ModelingToolkit array variable.
- `var`: symbolic array ModelingToolkit variable 
- `index`: index number 
```
"""
function to_latex_string(var)
    strvar, strindex = split_var_and_real_index(var)
    strvar = replace(strvar, "α" => "\\alpha")
    L"%$(strvar)_{%$(strindex)}(t)"
end  


"""
Replace ticks from x or y axis of a simple Plots with a LaTeXString version.
# Arguments 
- `a_plot::Plots.SubPlot`: a simple Plot, this not intended to work with plots with subplots.
- `axis::Symbol`: options are `:x` and `:y`. 
"""
function latexify_ticks!(a_plot, axis)
    if axis == :x  
        oldticks = Plots.xticks(a_plot)
    elseif axis == :y
        oldticks = Plots.yticks(a_plot)
    end 
    newticks = latexstring.(replace.(oldticks[2], "×" => "\\times"))
    if axis == :x  
        plot!(a_plot, xticks = (oldticks[1], newticks))
    elseif axis == :y
        plot!(a_plot, yticks = (oldticks[1], newticks))
    end 
end 

"""
Applies `latexify_ticks!` in x and y axis to every subplot of `a_plot`
"""
function latexify_ticks!(a_plot)
    for subplot in 1:length(a_plot)
        latexify_ticks!(a_plot[subplot], :y)
        if ! remove_xticks(subplot)
            latexify_ticks!(a_plot[subplot], :x)
        end
    end 
end

#===========================
Plot with scientific notation:
Conservé este código porque es útil para graficar el caso sintético.
===========================# 
#=
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
        annotate!(a_plot, (-0.07, 1.1), text(L"\times 10^%$(exponent)", :left, 8), subplot = subplot)
    end
end

vars_in_state(state) = ((state-1)*n + 1):((state-1)*n + n)
calculate_exponent(ymax) = floor(Int, log10(ymax)) 

function calculate_scaling_exponents(xs, state)
    ymax = maximum(xs[:,vars_in_state(state)])
    calculate_exponent(ymax)
end

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
    label = to_latex_string(symstates[index]) # antes era "s$index"
    plot!(ts, xs[:,index] * scaling_factor, label = label, ribbon = sqrt.(Ps[index,index,:]) * scaling_factor; kwargs...)
end 

remove_xticks(state) = ! (state in 5:6)

remove_xticks!(a_plot, subplot) = plot!(a_plot[subplot], xticks = (Plots.xticks(a_plot[subplot])[1], ["" ]))

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
function plot_all_states_grid(ts, xs, Ps, symstates)
    scaling_exponents = [calculate_scaling_exponents(xs, state) for state in 1:5]
    scaling_factors = 10 .^(-Float64.(scaling_exponents))
    a_plot = plot(layout=(3,2),framestyle=:box, link = :x, size = (800, 450));
    for state in 1:5 # estados 
        for i = 1:n # clases  
            plot_smoothed!(a_plot, ts, xs, Ps, symstates, (state-1)*n + i, scaling_factor = scaling_factors[state], subplot = state, fillalpha = 0.1) # 10^5 
        end 
        add_scix10_to_plot!(a_plot, state, scaling_exponents[state])
        if remove_xticks(state)
            remove_xticks!(a_plot, state)
        end
    end 
    plot_smoothed!(a_plot, ts, xs, Ps, symstates, 5*n + 1, subplot = 6, fillalpha = 0.1); # 10^5 
    if ! one_control 
        for i = 2:n # comunas 
            plot_smoothed!(a_plot, ts, xs, Ps, symstates, 5*n + i, subplot = 6, fillalpha = 0.1) # 10^5 
        end 
    end 
    a_plot 
end 