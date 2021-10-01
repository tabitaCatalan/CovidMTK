#=
Inner aux functions of model 
=#

using ModelingToolkit

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
    contact_rate = collect(TRM * collect(β .* (TRM' * E) ./ (TRM' * N)))
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

zero_dynamics(α) = [D(αi) ~ 0. for ai in α] 
function common_control(α, control_pieces, t) 
    [αi ~ control_pieces(t) for ai in α]
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

function get_alpha(one_control, t, n) 
    if one_control 
        @variables α(t) 
    else 
        @variables α[1:n](t)
    end 
    α
end 


zero_dynamics(α) = [D(αi) ~ 0. for αi in α] 

"""
# Arguments 
- `α`: Num, of Symbolics, dependent of `t` `α(t)`.
- `control_pieces`: registered in Symbolics function of t 
- `t`: Num, of Symbolics. Independent variable 
"""
function common_control(α, control_pieces, t) 
    [αi ~ control_pieces(t) for αi in α]
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
function epi_model_known_input(t, n, m, residence_times_data, control_pieces, one_control = true; name)
    α = get_alpha(one_control, t, n) 
    control_eqs = common_control(α, control_pieces, t) 
    eqs = eqs_common_epi_model(n, m, t, residence_times_data, α, control_eqs)
    ODESystem(eqs;name)
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
"""
function epi_model_unknown_input(t, n, m, residence_times_data, one_control = true; name)
    α = get_alpha(one_control, t, n)
    eqs = eqs_common_epi_model(n, m, t, residence_times_data, α, zero_dynamics(α))
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

post_process(str) = replace(replace(str, " " => ""), "." => "-")
function make_img_name(p)
    gamma_e_index = 1; gamma_i_index = 2; beta_2_index = 6; 
    gamma_e = last(p[gamma_e_index]); gamma_i = last(p[gamma_i_index]); beta_2 = last(p[beta_2_index])
    post_process("gamma_e_$gamma_e _gamma_i_$gamma_i _beta_2_$beta_2")
end
