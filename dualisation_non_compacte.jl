using JuMP
using CPLEX
#using GLPK



function robuste_dualisation_non_compacte(file_name::String)
    file_acces = "./data/" * file_name
    include(file_acces)


    # Taille de la grille
    n = size(t, 1)
    
    # Créer le modèle
    m = Model(CPLEX.Optimizer)

    #m = Model(GLPK.Optimizer)
    
       

    @variable(m, x[1:n, 1:n], Bin);
    @variable(m, 0 <= u[1:n]);

    @variable(m, lambda_1 >= 0 );
    @variable(m, lambda_2 >= 0 );
    @variable(m, alpha[1:n, 1:n]>= 0);
    @variable(m, beta[1:n, 1:n]>= 0);



    @objective(m,Min,sum(beta[i,j]*2 + alpha[i,j] for i in 1:n, j in 1:n) 
                     + sum(x[i,j]*t[i,j] for i in 1:n, j in 1:n)
                        + T*lambda_1 + T*T * lambda_2)




    @constraint(m, [i in 2:n] , sum(x[i,j] for j in 1:n ) == 1 )
    @constraint(m, [j in 2:n] , sum(x[i,j] for i in 1:n ) == 1 )
    @constraint(m, sum(x[1,j] for j in 2:n) == sum(x[i,1] for i in 2:n))
                    
                   
    @constraint(m, sum(x[i,i] for i in 1:n) == 0) # on ne peut pas aller de i à i

    
    @constraint(m, [i in 2:n,j in 2:n, i!=j], u[j] - u[i] >= d[i] - C*(1-x[i,j]))

   

    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] + th[j]) * x[i,j] <= lambda_1 + alpha[i,j] )
    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] * th[j]) * x[i,j] <= lambda_2 + beta[i,j] )


        # Ajout dynamique des contraintes de sous-tour (SEC)
    function add_subtour_constraints(m)
        sol = value.(x)  # Solution courante
        visited_nodes = [i for i in 2:n if sum(sol[i, j] for j in 1:n) > 0.5]

        if length(visited_nodes) > 1
            # Si sous-tour détecté, on l'interdit
            @constraint(m, sum(x[i, j] for i in visited_nodes, j in visited_nodes) <= length(visited_nodes) - 1)
        end
    end

    # Résolution avec ajout dynamique des SEC
    iteration = 0
    while true
        iteration += 1
        optimize!(m)
        
        prev_obj = objective_value(m)
        add_subtour_constraints(m)  # Ajout des SEC

        optimize!(m)  
        new_obj = objective_value(m)

        if abs(prev_obj - new_obj) < 1e-5  # Si convergence
            break
        end
    end


    optimize!(m)
    return JuMP.objective_value(m)

end
"""
for tour in sous_tour
        complément = setdiff(1:n, tour)
        @constraint(m, sum(x[i,j] for i in tour, j in complément) >= sum(d[i] for i in tour) / C)
    end
"""

list_name = ["n_5-euclidean_false","n_5-euclidean_true",
                "n_6-euclidean_false","n_6-euclidean_true",
                "n_7-euclidean_false","n_7-euclidean_true",
                "n_8-euclidean_false","n_8-euclidean_true",
                "n_9-euclidean_false","n_9-euclidean_true",
                "n_10-euclidean_false","n_10-euclidean_true",
                "n_11-euclidean_false","n_11-euclidean_true",
                "n_12-euclidean_false","n_12-euclidean_true",
                "n_13-euclidean_false","n_13-euclidean_true",
                "n_14-euclidean_false","n_14-euclidean_true",
                "n_15-euclidean_false","n_15-euclidean_true",
                "n_16-euclidean_false","n_16-euclidean_true",]
liste_reduite = ["n_5-euclidean_true"]

function main(list_name)
    results = []
    global_start = time()
    for file_name in list_name
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
  