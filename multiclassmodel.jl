#=
Inner aux functions of model 
=#

using ModelingToolkit

#=============================================
Generic functions 
=============================================#

"""
- `γₑ::Num` symbolic variable of MTK, time dependent 
- `gamma_e::Float64` a constant value for `γₑ`
"""
function constant_value_equation(sym, value, t) 
    [sym ~ value]
end 

zero_dynamics(α, t) = [Differential(t)(αi) ~ 0. for αi in α] 


#===============================================# 
function common_dynamics(D, x, p, λ)
    γₑ, γᵢ = p 
    S, E, I, R, C, N = x  
    n = length(S) 
    [
        [D(S[i]) ~ - λ[i] .* S[i] for i in 1:n];
        [D(E[i]) ~ λ[i] .* S[i] - γₑ * E[i] for i in 1:n];
        [D(I[i]) ~ γₑ * E[i] - γᵢ * I[i] for i in 1:n]; 
        [D(R[i]) ~ γᵢ * I[i] for i in 1:n]; 
        [D(C[i]) ~ γₑ * E[i] for i in 1:n];
    ]
end 

function environment_contact_rate(t, x, λ, TRM, residence_times_data, α, β)
    @variables i j 

    n,m = size(TRM) 
    S, E, I, R, C, N = x 
    contact_rate = collect(TRM * collect(β .* (TRM' * I) ./ (TRM' * N)))
    if length(α) == 1 
        contact_rate = λ .~ α * contact_rate
    else 
        contact_rate = λ .~ (α .* contact_rate)
    end 
    [
        vec([TRM[i,j] ~ residence_times_data(t, i, j) for i in 1:n, j in 1:m]); # t + 28 hay que agregarlo en algún lado 
        collect(contact_rate);
    ]
end 


function eqs_common_epi_model(n, m, t, residence_times_data, α, control_eqs)
    a = 1
    @parameters γₑ γᵢ β[1:m] N[1:n]
    @variables S[1:n](t) E[1:n](t) I[1:n](t) R[1:n](t) C[1:n](t) λ[1:n](t)
    @variables TRM[1:n, 1:m](t)
     

    D = Differential(t) 

    #S, E, I, R, C, N = x 
    x = (S, E, I, R, C, N) 

    [
        environment_contact_rate(t, x, λ, TRM, residence_times_data, α, β);
        common_dynamics(D, x, (γₑ, γᵢ), λ);
        control_eqs;
    ];
end 

function eqs_epi_model_variable_rates(n, m, t, residence_times_data, α, control_eqs, (γₑ, γᵢ), rate_eqs)
    a = 1
    @parameters β[1:m] N[1:n]
    @variables S[1:n](t) E[1:n](t) I[1:n](t) R[1:n](t) C[1:n](t) λ[1:n](t)
    @variables TRM[1:n, 1:m](t)
     

    D = Differential(t) 

    #S, E, I, R, C, N = x 
    x = (S, E, I, R, C, N) 

    [
        environment_contact_rate(t, x, λ, TRM, residence_times_data, α, β);
        common_dynamics(D, x, (γₑ, γᵢ), λ);
        control_eqs;
        rate_eqs;
        # zero_dynamics(γₑ); para el caso unknown 
    ];
end 

function get_alpha(one_control, t, n) 
    if one_control 
        @variables α(t) 
    else 
        @variables α[1:n](t)
    end 
    α
end 

function get_rates(variable_rate, t)
    if variable_rate
        @variables γₑ(t) γᵢ(t)
    else
        @parameters γₑ γᵢ
    end  
    (γₑ, γᵢ)
end 

function get_rate_equations(variable_rate::Bool, (γₑ, γᵢ), gamma_e, gamma_i, t; known::Bool)
    if variable_rate
        if known
            rate_eqs = [constant_value_equation(γₑ, gamma_e, t); constant_value_equation(γᵢ, gamma_i, t)]
        else
            rate_eqs = zero_dynamics([γₑ, γᵢ], t)
        end
    else
        rate_eqs = []
    end  
    rate_eqs
end 



"""
# Arguments 
- `α`: Num, of Symbolics, dependent of `t` `α(t)`.
- `control_pieces`: registered in Symbolics function of t 
- `t`: Num, of Symbolics. Independent variable 
"""
function common_control(α, control_pieces, t) 
    [αi ~ control_pieces(t) for αi in α]
end 

function two_control(α, control, t) 
    [α[i] ~ control(t,i) for i in 1:2]
end 


#=
Principal model functions 
=#

"""
# Arguments 
- `t::Num`: Symbolics independent variable.
- `residence_times_data`: registered in Symbolics function of `(t, i, j)`.
- `control_pieces`: registered in Symbolics function of `t`. 
- `one_control::Bool`. If `true` a common control `α(t)` is used for each class. If `false`, 
    a different control for each class `αᵢ(t)` is used instead.
- `t`: Num, of Symbolics. Independent variable 
"""
function epi_model_known_input(t, n, m, residence_times_data, control_pieces, one_control = true, variable_rate = false, gamma_e = 1/5, gamma_i = 1/7; name)
    α = get_alpha(one_control, t, n) 
    if one_control
        control_eqs = common_control(α, control_pieces, t) 
    else 
        control_eqs = two_control(α, control_pieces, t) 
    end

    γₑ, γᵢ = get_rates(variable_rate, t)
    rate_eqs = get_rate_equations(variable_rate, (γₑ, γᵢ), gamma_e, gamma_i, t, known = true)
    eqs = eqs_epi_model_variable_rates(n, m, t, residence_times_data, α, control_eqs, (γₑ, γᵢ), rate_eqs)
    ODESystem(eqs;name)
end 

α = get_alpha(one_control, t, n) 
if one_control
    control_eqs = common_control(α, control_pieces, t) 
else 
    control_eqs = two_control(α, controls, t) 
end

"""
# Arguments 
- `t::Num`: Symbolics independent variable.
- `n`: number of classes 
- `m`: number of environments 
- `residence_times_data`: registered in Symbolics function of ``(t, i, j)`` where ``t`` is time and ``(i,j)`` are
    the indexes of a matrix, ``1 \\leq i \\leq n, 1 \\leq j \\leq m``.
- `one_control::Bool`. If `true` a common control `α(t)` is used for each class. If `false`, 
    a different control for each class `αᵢ(t)` is used instead.
- `variable_rate::Bool`: if `true`, then `γₑ` and `γᵢ` are considered as part of the augmented state.
- `gamma_e = 1/5`: only used if `variable_rate` is `true`.
"""
function epi_model_unknown_input(t, n, m, residence_times_data, one_control = true, variable_rate = false, gamma_e = 1/5, gamma_i = 1/7; name)
    α = get_alpha(one_control, t, n)
    γₑ, γᵢ = get_rates(variable_rate, t)
    rate_eqs = get_rate_equations(variable_rate, (γₑ, γᵢ), gamma_e, gamma_i, t, known = false)
    eqs = eqs_epi_model_variable_rates(n, m, t, residence_times_data, α, zero_dynamics(α, t), (γₑ, γᵢ), rate_eqs)
    ODESystem(eqs;name)
end 

#=
Aux functions to create initial conditions and parameters 
=#
"""
Devuelve una tupla `(n,m)`, donde  `n` es el número de clases que usa el sistema 
y `m` el número de ambientes.
# Example 
```
julia> get_size(episys_known)
(2,2)
```
"""
get_size(system) = size(system.TRM)

function make_x_k(system, S0, E0, I0, R0, C0)
    n, m = get_size(system)
    [
        [S[i] => S0[i] for i in 1:n]; 
        [E[i] => E0[i] for i in 1:n];
        [I[i] => I0[i] for i in 1:n]; 
        [R[i] => R0[i] for i in 1:n];
        [C[i] => C0[i] for i in 1:n];
    ]    
end 

function make_alpha(system, α0) 
    n, m = get_size(system)
    
    if length(system.α) > 1
        α = get_alpha(false, t, n) 
        if length(α0) == n 
            control_eq = [α[i] => α0[i] for i in 1:n]
        elseif length(α0) == 1
            control_eq = [α[i] => α0 for i in 1:n]
        end 
    else 
        α = get_alpha(true, t, n) 
        control_eq = [α => α0] 
    end 
    control_eq
end 

function make_rate(gamma_e, gamma_i)
    [γₑ => gamma_e, γᵢ => gamma_i]
end 

function make_x_uk(system, α0, S, E, I, R, C)
    
    control_eq = make_alpha(system, α0) 
    [
        control_eq;
        make_x_k(system, S, E, I, R, C)
    ]    
end 

function make_p(system, gammae, gammai, total, beta)
    n, m = get_size(system)
    [
        [γₑ => gammae, γᵢ => gammai]; # midiendo en semanas... un lío con la matrix P 
        [N[i] => total[i] for i = 1:n];
        [β[i] => beta[i] for i = 1:m];
    ]
end 



#=
File names generation 
=#
using Printf: @sprintf 
post_process(str) = replace(replace(str, " " => ""), "." => "-")
function make_img_name(p)
    gamma_e_index = 1; gamma_i_index = 2; beta_2_index = 6; 
    gamma_e = last(p[gamma_e_index]); gamma_i = last(p[gamma_i_index]); beta_2 = last(p[beta_2_index])
    str = @sprintf "gamma_e_%.2f _gamma_i_%.2f _beta_2_%.2f" gamma_e gamma_i beta_2 
    post_process(str)
end
