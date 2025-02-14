using JuMP
using CPLEX
#using GLPK



using JuMP, CPLEX

function cvrp_non_compact(file_name::String)
    # Chargement des données
    file_acces = "./data/" * file_name
    include(file_acces)  # Fichier doit contenir t, d, C

    n = size(t, 1)  # Nombre de clients + 1 (dépôt)

    # Modèle JuMP
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", 180)

    # Variables de décision
    @variable(m, x[1:n, 1:n], Bin)   # x[i,j] = 1 si l’arc (i,j) est utilisé
    @variable(m, 0 <= y[1:n, 1:n])   # y[i,j] = charge transportée sur l’arc (i,j)

    # Objectif : Minimiser la distance totale parcourue
    @objective(m, Min, sum(t[i,j] * x[i,j] for i in 1:n, j in 1:n))

    # Chaque client est visité exactement une fois
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, [j in 2:n], sum(x[i,j] for i in 1:n if i != j) == 1)

    # Contrainte de flot (capacité des véhicules)
    @constraint(m, [i in 2:n, j in 2:n,  i != j],  y[i,j] >= d[j] * x[i,j])  # Si un arc est pris, il doit transporter au moins la demande
    @constraint(m, [i in 2:n, j in 2:n,  i != j],  y[i,j] <= (C - d[i]) * x[i,j])  # Le véhicule ne doit pas dépasser sa capacité

    # Conservation du flot
    @constraint(m, [i in 2:n], 
        sum(y[j,i] for j in 1:n if j != i) == d[i] + sum(y[i,j] for j in 1:n if j != i))

    # Les véhicules partent et reviennent au dépôt
    @constraint(m, sum(x[1,j] for j in 2:n) == sum(x[i,1] for i in 2:n))

    # Suppression des auto-boucles
    @constraint(m, [i in 1:n], x[i,i] == 0)

    # Ajout dynamique des contraintes de sous-tour
    function add_subtour_constraints()
        sol_x = value.(x)  
        visited_nodes = [i for i in 2:n if sum(sol_x[i, j] for j in 1:n) > 0.5]

        if length(visited_nodes) > 1
            @constraint(m, sum(x[i, j] for i in visited_nodes, j in visited_nodes) <= length(visited_nodes) - 1)
        end
    end

    # Boucle d’optimisation avec suppression des sous-tours
    iteration = 0
    while true
        iteration += 1
        optimize!(m)
        
        prev_obj = objective_value(m)
        add_subtour_constraints()  # Ajout des SEC

        optimize!(m)  
        new_obj = objective_value(m)

        if abs(prev_obj - new_obj) < 1e-5  # Convergence atteinte
            break
        end
    end

    return objective_value(m)
end



list_name = ["n_5-euclidean_false","n_5-euclidean_true",
                "n_6-euclidean_false","n_6-euclidean_true",
                "n_7-euclidean_false","n_7-euclidean_true",
                "n_8-euclidean_false","n_8-euclidean_true",
                "n_9-euclidean_false","n_9-euclidean_true",
                "n_10-euclidean_false","n_10-euclidean_true",
                "n_11-euclidean_false","n_11-euclidean_true",
                "n_12-euclidean_false","n_12-euclidean_true",
                "n_13-euclidean_false","n_13-euclidean_true",
                "n_14-euclidean_false"
                ]
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
    fout = open("output_non_compact.txt", "w")
    # Ecrire "test" dans ce fichier
    println(fout, results)
    close(fout)
end
  