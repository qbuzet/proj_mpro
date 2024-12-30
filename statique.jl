using JuMP
using CPLEX
using GLPK


function statique(file_name::String)
    file_acces = "./data/" * file_name
    include(file_acces)

    # Créer le modèle
    m = Model(CPLEX.Optimizer)
    #m = Model(GLPK.Optimizer)

    @variable(m, x[1:n, 1:n], Bin);
    @variable(m, 0 <= u[1:n] );

    @objective(m,Min, sum(x[i,j]*t[i,j] for i in 1:n, j in 1:n))



    @constraint(m, [i in 2:n] , sum(x[i,j] for j in 1:n ) == 1 )
    @constraint(m, [j in 2:n] , sum(x[i,j] for i in 1:n ) == 1 )
    @constraint(m, sum(x[1,j] for j in 2:n) == sum(x[i,1] for i in 2:n))

    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n], u[i] <= C*(1-x[1,i]))

    @constraint(m, [i in 2:n,j in 2:n, i!=j], u[j] - u[i] >= d[i] - C*(1-x[i,j]))
    @constraint(m, sum(x[i,i] for i in 1:n) == 0) # on ne peut pas aller de i à i
    

    optimize!(m)
    println("Valeur de l’objectif : ", JuMP.objective_value(m))

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        
end
end