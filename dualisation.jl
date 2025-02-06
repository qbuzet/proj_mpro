using JuMP
#using CPLEX
using GLPK

function robuste_dualisation(file_name::String)
    file_acces = "./data/" * file_name
    include(file_acces)

    # Taille de la grille
    n = size(t, 1)

    # Créer le modèle
    #m = Model(CPLEX.Optimizer)
    m = Model(GLPK.Optimizer)

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
                    
    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n], u[i] <= C*(1-x[1,i]))
                   
    @constraint(m, [i in 2:n,j in 2:n, i!=j], u[j] - u[i] >= d[i] - C*(1-x[i,j]))
    @constraint(m, sum(x[i,i] for i in 1:n) == 0) # on ne peut pas aller de i à i
                        
   

    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] + th[j]) * x[i,j] <= lambda_1 + alpha[i,j] )
    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] * th[j]) * x[i,j] <= lambda_2 + beta[i,j] )

    optimize!(m)
    return JuMP.objective_value(m)

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
                "n_14-euclidean_false","n_14-euclidean_true",
                "n_15-euclidean_false","n_15-euclidean_true",
                "n_16-euclidean_false","n_16-euclidean_true",
                "n_17-euclidean_false","n_17-euclidean_true",
                "n_18-euclidean_false","n_18-euclidean_true",
                "n_19-euclidean_false","n_19-euclidean_true",
                "n_20-euclidean_false","n_20-euclidean_true",
                "n_25-euclidean_false","n_25-euclidean_true",
                "n_30-euclidean_false","n_30-euclidean_true",
                "n_35-euclidean_false","n_35-euclidean_true",
                "n_40-euclidean_false","n_40-euclidean_true",
                "n_45-euclidean_false","n_45-euclidean_true",
                "n_50-euclidean_false","n_50-euclidean_true",
                "n_55-euclidean_false","n_55-euclidean_true",
                "n_60-euclidean_false","n_60-euclidean_true",
                "n_65-euclidean_false","n_65-euclidean_true",
                "n_70-euclidean_false","n_70-euclidean_true",
                "n_75-euclidean_false","n_75-euclidean_true",
                "n_80-euclidean_false","n_80-euclidean_true",
                "n_85-euclidean_false","n_85-euclidean_true",
                "n_90-euclidean_false","n_90-euclidean_true",
                "n_95-euclidean_false","n_95-euclidean_true",
                "n_100-euclidean_false","n_100-euclidean_true",
                ]
liste_reduite = ["n_5-euclidean_false"]
function main(liste)
    results = []
    global_start = time()
    for file_name in liste
        if occursin("euclidean", file_name)
            start = time()           
            println("Fichier : ", file_name)
            obj = robuste_dualisation(file_name)
            exec_time = time() - start
            global_exec_time = time() - global_start
            push!(results, (file_name, obj, exec_time, global_exec_time))
            if global_exec_time > 3600
                break
            end
        end
    end
    fout = open("output_dualisation.txt", "w")
    # Ecrire "test" dans ce fichier
    println(fout, results)
    close(fout)
end
  