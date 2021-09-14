#=
Linear coefficients 
=#

#= Assume we jave an ODESystem `epi_system`. We can get 
an array of its states by using `states(epy_system)`.
We have also a vector of expressiones 
=#

function get_linear_coefficients(sym_expression, var, system)
    sysstates = states(system)
    identifmap = [state => isequal(var, state) ? 1. : 0. for state in sysstates]
    vec = ModelingToolkit.varmap_to_vars(identifmap, sysstates)
     # expression=Val{false} permite generar la función durante RunTime 
     # sin eso da un error raro, relacionado como con los scopes de las cosas
    eval(build_function(sym_expression, sysstates,  expression=Val{false}))(vec)
end 

get_linear_coefficients(2C[1] + E[1], C[1], simple_epi_system) == 2.

function get_observacion_map(sym_expression, system)
    [state => get_linear_coefficients(sym_expression, state, system) for state in states(system)]
end 


isequal(get_observacion_map(C[1] + I[2], simple_epi_system), [
    S[1] => 0., 
    S[2] => 0., 
    E[1] => 0., 
    E[2] => 0., 
    R[1] => 0.,
    R[2] => 0., 
    I[1] => 0.,
    I[2] => 1.,
    C[1] => 1.,
    C[2] => 0.,
    α[1] => 0., 
    α[2] => 0.]) # da false porque no compara bien...



get_observation_vector(sym_exp, system) = ModelingToolkit.varmap_to_vars(get_observacion_map(sym_exp, system), states(system))'

get_observation_vector(C[1]+2R[2], simple_epi_system)

get_observacion_matrix(exp_vec, system) = vcat([get_observation_vector(exp, system) for exp in exp_vec]...)
get_observacion_matrix([C[1], C[2]], simple_epi_system) 
