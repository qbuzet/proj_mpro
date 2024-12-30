using JuMP
#using CPLEX
using GLPK


function statique(file_name::String)
    file_acces = "./data/" * file_name
    include(file_acces)

    # Créer le modèle
    #m = Model(CPLEX.Optimizer)
    m = Model(GLPK.Optimizer)

    @variable(m, x[1:n, 1:n, 1:n], Bin);
    @variable(m, y[1:n], Bin);
    @variable(m, 0 <= u[1:n, 1:n] <= n);

    @objective(m,Min, sum(x[i,j,k]*t[i,j] for i in 1:n, j in 1:n, k in 1:n))


    @constraint(m, [k in 1:n] , sum(x[i,j,k] * d[i] for i in 1:n, j in 1:n ) <= C*y[k] )

    @constraint(m, [i in 2:n] , sum(x[i,j,k] for k in 1:n, j in 1:n ) == 1 )
    @constraint(m, [j in 2:n] , sum(x[i,j,k] for k in 1:n, i in 1:n ) == 1 )

    @constraint(m, [i in 1:n,j in 1:n,k in 1:n], u[j,k]>= u[i,k] + 1 - n*(1-x[i,j,k]))

    @constraint(m, sum(x[1,j,k] for k in 1:n, j in 2:n) == sum(y[k] for k in 1:n))
    @constraint(m, sum(x[i,1,k] for k in 1:n, i in 2:n) == sum(y[k] for k in 1:n))
    @constraint(m,[k in 1:n-1], y[k] >= y[k+1])
    @constraint(m, sum(x[1,1,k] for k in 1:n) == 0)

    optimize!(m)

    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        println("Valeur de l’objectif : ", JuMP.objective_value(m))
        println("Nombre de véhicules : ", JuMP.value(sum(y[k] for k in 1:n)))
        # Reconstruction des chemins
       for k in 1:n
            i = 1
            j0 = 0
            while j0 != 1
                for j in 1:n
                    if value(x[i,j,k]) == 1
                        println("Véhicule ", k, " : ", i, " -> ", j)
                        i = j
                        j0 = j
                        break
                    end
                
                
            end
        end
    end
end
end