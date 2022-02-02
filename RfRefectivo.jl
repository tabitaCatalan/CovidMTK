#=
Cálculo del R efectivo 
=#

# episys_known es el sistema que me interesa 
# Obtener ecuaciones del sistema 
# Obtener entradas y salidas de compartimientos sanos 

generate_function(simple_episys_known) 


γᵢ .* I 



in_flow =     [
    [0 for i in 1:n]; 
    [γᵢ * I[i] for i in 1:n]; 
]

Symbolics.jacobian(in_flow ,[S;R]) 

out_flow = [
    [λ[i] .* S[i] for i in 1:n];
    [0 for i in 1:n]; 
]

eqs_in, eqs_out = epi_model_known_flows(t, n, m, adjusted_rtm, control_pieces, false, variable_rate,gamma_e_real, gamma_i_real, variable_beta, beta_exterior_real)


@named in_flow_sys = ODESystem(eqs_in) 
in_flow_sys = structural_simplify(in_flow_sys)

subs = Dict([λ[1] => eqs_in[11].rhs])
subs = Dict(λ[1] => I[1])
substitute(out_flow, subs)

Symbolics.jacobian([eqs_in[11].rhs] ,[I;R]) 
@variables xvar yvar

substitute(xvar^2 - yvar, Dict(xvar => yvar))

@nonamespace episys_known.α 


calculate_jacobian(in_flow_sys, sparse = true) # esto parece que sirve ... 

jac2 = generate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
#jac = generate_jacobian(simple_episys_uknown); # esto debería ser la función que devuelve el jacobiano 

#jac2 = generate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
system_jacobian = eval(jac2[1]);

#jac2 = calculate_jacobian(simple_episys_uknown, sparse = true); # esto debería ser la función que devuelve el jacobiano 
#jac1 = calculate_jacobian(simple_episys_uknown);
pvec = ModelingToolkit.varmap_to_vars(p,parameters(simple_episys_uknown));
u0vec = ModelingToolkit.varmap_to_vars(u0,states(simple_episys_uknown));



function flow_factory()
    @variables t S(t) E(t) I(t) R(t)
    x = [S, R]; y = [E, I] 
    D = Differential(t) 
    [D(x) .~ 0.]
end 

#===================Composable Model for Rt =======================#

@variables t 

"""
- `residence_times_data`: function of `(t, i, j)` con `i` clase y `j` ambiente.
"""
function TimeResidenceMatrix(residence_times_data, n, m; name)
    @variables TRM[1:n, 1:m](t)
    eqs = vec([TRM[i,j] ~ residence_times_data(t, i, j) for i in 1:n, j in 1:m]); 
    ODESystem(eqs, t; name=name)
end 


ODESystem(collect(D.(S) .~ 0), t; name = :S0)

function MultiClassSanitaryFactor(control_function; name)
    @variables α[1:n](t)
    control_eqs = [α[i] ~ control_function(t, i) for i in 1:n]
    ODESystem(control_eqs, t; name=name)
end

function EnvironmentalRisk(variable_beta, beta_in, beta_ext; name) 
    @parameters βᵢₙ = beta_in 
    if variable_beta
        @variables βₑₓ(t) = beta_ext
        ODESystem(Equation[], t, [βₑₓ], [βᵢₙ]; name=name)
    else 
        @parameters βₑₓ = beta_ext
        ODESystem(Equation[], t, [], [βₑₓ, βᵢₙ]; name=name)
    end 
end

@variables ξ(t) = 26.
EnvironmentalRisk(false, 1., beta_exterior; name = :beta)

function Susceptible(n; name)
    @variables S[1:n](t)
    ODESystem(Equation[], t, S, []; name=name)
end 

function Exposed(n; name)
    @variables E[1:n](t)
    ODESystem(Equation[], t, E, []; name=name)
end 


function Infected(n; name)
    @variables I[1:n](t)
    ODESystem(Equation[], t, I, []; name=name)
end 


function Recovered(n; name)
    @variables R[1:n](t)
    ODESystem(Equation[], t, R, []; name=name)
end 

function ContactRate(residence_times_data, control_function, n, m, beta_ext; name)
    @variables λ[1:n](t) 
    @named RR = TimeResidenceMatrix(residence_times_data, n, m)
    @named alpha = MultiClassSanitaryFactor(control_function)
    @named beta = EnvironmentalRisk(true, 1., beta_ext)
    @named infected = Infected(n) 
    @named susc = Susceptible(n) 
    @parameters N[1:n]
    β = [beta.βᵢₙ, beta.βₑₓ]
    contact_rate = collect(RR.TRM * collect(β .* (RR.TRM' * infected.I) ./ (RR.TRM' * N)))
    lambda_eqs = collect.(λ .~ alpha.α .* contact_rate .* susc.S)
    compose(ODESystem(lambda_eqs, t; name=name), RR, alpha, beta, infected, susc)
end 

ContactRate(adjusted_rtm, comportamiento_normal, n, m, beta_exterior; name = :my_cr)

function InFlux(n; name)
    @variables IS[1:n](t) IR[1:n](t) 
    @parameters γᵢ 

    @named infected = Infected(n) 

    D = Differential(t)
    in_eqs = [
        [D(IS[i]) ~ 0 for i in 1:n]; 
        [D(IR[i]) ~ γᵢ * infected.I[i] for i in 1:n]; 
    ]
    compose(ODESystem(in_eqs, t; name=name), infected)
end 
function OutFlux(residence_times_data, control_function, n, m, beta_ext; name)
    @variables OS[1:n](t) OR[1:n](t) 
    D = Differential(t)

    @named lambda = ContactRate(residence_times_data, control_function, n, m, beta_ext)
    out_eqs = [
        [D(OS[i]) ~ lambda.λ[i] for i in 1:n];
        [D(OR[i]) ~ 0 for i in 1:n]; 
    ]
    compose(ODESystem(out_eqs, t; name=name), lambda)
end 

function DiseaseFreeDynamics(residence_times_data, control_function, n, m, beta_ext; name)
    @named out_flux = OutFlux(residence_times_data, control_function, n, m, beta_ext)
    @named in_flux = InFlux(n)
    
    @named susc = Susceptible(n) 
    @named rec = Recovered(n) 
    D = Differential(t)
    eqs = [
        collect(D.(susc.S) .~ in_flux.IS - out_flux.OS);
        collect(D.(rec.R) .~ in_flux.IR - out_flux.OR);
        # connect susc.S con out_flux.lambda.susc.S 
        collect(susc.S .~ out_flux.lambda.susc.S);
    ]
    compose(ODESystem(eqs, t; name=name), out_flux, in_flux, susc, rec)
end 

function DiseaseDynamics(residence_times_data, control_function, n, m, beta_ext; name)
    
    @named expo = Exposed(n) 
    @named inf = Infected(n) 
    D = Differential(t)

    @named lambda = ContactRate(residence_times_data, control_function, n, m, beta_ext)

    @parameters γᵢ γₑ
    
    eqs = [
        [D(expo.E[i]) ~ lambda.λ[i] - γₑ * expo.E[i] for i in 1:n];
        [D(inf.I[i]) ~ γₑ * expo.E[i] - γᵢ * inf.I[i] for i in 1:n];
        # connect lambda.infected.I con inf.I 
        #connect(lambda.infected.I, inf.I);
        collect(lambda.infected.I .~ inf.I);
    ]
    compose(ODESystem(eqs, t; name=name), lambda, expo, inf)
end 

# residence_times_data = adjusted_rtm
# control_function = comportamiento_normal  


@named disease_dyn = DiseaseDynamics(adjusted_rtm, comportamiento_normal, n, m, beta_exterior)
@named diseasefree_dyn = DiseaseFreeDynamics(adjusted_rtm, comportamiento_normal, n, m, beta_exterior)


