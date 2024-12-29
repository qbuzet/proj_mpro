using JuMP
#using CPLEX
using GLPK

function robuste_dualisation(file_name::String)
    file_acces = "./data/" + file_name
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

    @constraint(m, [i in 1:n,j in 1:n,k in 1:n], u[j,k]>= u[i,k] + 1 - n*(1-x[i,j,k]))

    @constraint(m, sum(x[1,j,k] for k in 1:n, j in 2:n) == sum(y[k] for k in 1:n))
    @constraint(m, sum(x[i,1,k] for k in 1:n, i in 2:n) == sum(y[k] for k in 1:n))

    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] + th[j]) * sum(x[i,j,k] for k in 1:n) >= lambda_1 + alpha[i,j] )
    @constraint(m, [j in 1:n, i in 1:n, i!=j], (th[i] * th[j]) * sum(x[i,j,k] for k in 1:n) >= lambda_2 + beta[i,j] )


  