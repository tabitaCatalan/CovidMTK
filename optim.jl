#=
Optimización de parámetros usando BlackBoxOptim 
=#

using BlackBoxOptim: bboptimize, best_candidate, best_fitness

function loss_function(x; rango = 1:400)
    try 
        a0 = x[1]
        u0 = initial_u0(a0)
        p = create_p(x)
        results, xs, Ps = kalman_iteration(u0, p)
        
        println("")
        #loss(results.analysis, observaciones, rango)
        loss(xs, observaciones, rango)
    catch 
        println("Error in $x")
        Inf 
    end
end 

#loss_function([0.011, 60., 0.19607843137254904, 0.125])
#loss_function([0.024567385731049017, 13.13089162772159,  0.10382086111103203, 0.31788134378524546])
#loss_function([0.016728510433320562, 34.87213978597166, 0.1451405896110646, 0.16229284534863325])

#loss(xs, observaciones, 1:400)
# a0 = x[1]; beta2 = x[2]; gamma_e = x[3]; gamma_i = x[4]

# --------> este es el que sirve... pero no lo quiero correr porque muere todo xD 
#opt = bboptimize(loss_function; SearchRange = [(1e-3, 3e-2), (1., 100.), (1/7, 1/3), (1/14, 1/5)]) 

#=

@inbounds x1, x2, x3, x4 = [0.024567385731049017,  13.13089162772159,   0.10382086111103203,   0.31788134378524546]

best_candidate(opt)
  0.024567385731049017
  13.13089162772159
  0.10382086111103203
  0.31788134378524546

best_fitness(opt)
113826.8206156808
=#
#=
best_candidate(opt)
4-element Array{Float64,1}:
  0.016728510433320562
 34.87213978597166
  0.1451405896110646
  0.16229284534863325

best_fitness(opt)
133239.68624681467
=#


#using Optim

