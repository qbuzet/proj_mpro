using JuMP
#using CPLEX
using GLPK

function robuste_dualisation_non_compacte(file_name::String)
        file_acces = "./data/" * file_name
        include(file_acces)


        # Taille de la grille
        n = size(t, 1)
    
        # Créer le modèle
        #m = Model(CPLEX.Optimizer)
        m = Model(GLPK.Optimizer)
    
        @variable(m, x[1:n, 1:n, 1:n], Bin);
        @variable(m, y[1:n], Bin);
        @variable(m, 0 <= u[1:n, 1:n] <= n, Int);
    
        @variable(m, lambda_1 >= 0 );
        @variable(m, lambda_2 >= 0 );
        @variable(m, alpha[1:n, 1:n]>= 0);
        @variable(m, beta[1:n, 1:n]>= 0);
    
    
    
        @objective(m,Min,sum(beta[i,j]*2 + alpha[i,j] for i in 1:n, j in 1:n) 
                         + sum(x[i,j,k]*t[i,j] for i in 1:n, j in 1:n, k in 1:n)
                            + T*lambda_1 + T*T * lambda_2)
    
    
    
        @constraint(m, [k in 1:n] , sum(x[i,j,k] * d[i] for i in 1:n, j in 1:n ) <= C*y[k] )
    
        @constraint(m, [i in 2:n] , sum(x[i,j,k] for k in 1:n, j in 1:n ) == 1 )
        @constraint(m, [j in 2:n] , sum(x[i,j,k] for k in 1:n, i in 1:n ) == 1 )
    
        @constraint(m, [i in 2:n,j in 2:n,k in 1:n, i!=j], u[j,k]>= u[i,k] + 1 - n*(1-x[i,j,k]))
        @constraint(m, [k in 1:n, j in 2:n], u[j,k]<=n*(1-x[1,j,k]))

        @constraint(m, sum(x[j,j,k] for k in 1:n, j in 1:n) == 0)

    
        @constraint(m, sum(x[1,j,k] for k in 1:n, j in 2:n) == sum(y[k] for k in 1:n))
        @constraint(m, sum(x[i,1,k] for k in 1:n, i in 2:n) == sum(y[k] for k in 1:n))
    
        @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] + th[j]) * sum(x[i,j,k] for k in 1:n) <= lambda_1 + alpha[i,j] )
        @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] * th[j]) * sum(x[i,j,k] for k in 1:n) <= lambda_2 + beta[i,j] )

        @constraint(m, [k in 1:n-1], y[k] >= y[k+1])

    optimize!(m)
    println("Valeur de l’objectif : ", JuMP.objective_value(m))

end
    


function main()
    results = []
    global_start = time()
    for file_name in ["n_6-euclidean_false"]
        if occursin("euclidean", file_name)
            start = time()           
            println("Fichier : ", file_name)
            obj = robuste_dualisation_non_compacte(file_name)
            exec_time = time() - start
            global_exec_time = time() - global_start
            push!(results, (file_name, obj, exec_time, global_exec_time))
        end
    end
    println(results)
end
  