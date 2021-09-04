#=
Cómo calcular el jacobiano de un sistema
=#
using ModelingToolkit, OrdinaryDiffEq 

@parameters t α γₑ
@variables S[1:2](t) 


simple_eqs = [
    D(S[1]) ~ -γₑ * S[1]^2
    D(S[2]) ~ - α * S[2]
]

sys_simple_eqs = ODESystem(simple_eqs, t) 

# esto da una expresión simbólica del jacobiano del lado derecho 
calculate_jacobian(sys_simple_eqs) 

#=
esto da un código de bajo nivel que puede usarse para construir 
una función que entregue el jacobiano 
`f_exp_jac` es una tupla. 
`f_exp_jac[1]` es una función de la forma
(vector de estados, vector de parametros, t(iempo)). 
`f_exp_jac[2]` es una función que muta el primer argumento,
de la forma 
(vector para allocar, vector de estados, vector de parametros, t(iempo)). 
=# 
f_exp_jac = generate_jacobian(sys_simple_eqs)

f_jac = eval(f_exp_jac[1])

x = [
    S[1] => 4.,
    S[2] => 2.
]
p₁ = [
    α => 0.01,
    γₑ => 1/5
] 
p1vec = ModelingToolkit.varmap_to_vars(p₁,parameters(sys_simple_eqs));
xvec = ModelingToolkit.varmap_to_vars(x,states(sys_simple_eqs));

f_jac(xvec, p1vec, 0.1 ) 


parameters(sys_simple_eqs)
dynamics = eval(generate_function(sys_simple_eqs)[1]) 
dynamics(xvec, p1vec, 0.)
