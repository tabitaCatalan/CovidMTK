#=
Aux functions for plotting 
=#
using Plots: annotate!, text, scatter, scatter! 
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
    longstrvar = string(var)
    leftbracketindex = findfirst('[', longstrvar)
    rightbracketindex = findfirst(']', longstrvar)
    strvar = longstrvar[1:leftbracketindex-1]
    strindex = longstrvar[leftbracketindex+1:rightbracketindex-1]
    strvar, strindex
end 
- `var`: symbolic array ModelingToolkit variable 
# Example 
```julia
@variables t S[1,2](t)
string(S[1]) #"S[1](t)"
to_string(S, 1) #"S₁(t)"
```
"""
function to_string(var, index)
    longstrvar = string(var[index])
    strvar = longstrvar[1:findfirst('[', longstrvar)-1]
    strvar * to_subindex(string(index)) * "(t)"
end  


#===========================
Plot with scientific notation 
===========================# 
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


