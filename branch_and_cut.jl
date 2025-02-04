using JuMP
using CPLEX

file = "data/instance_n5.txt"
include(file)

function solve_SP(xx)
    sm = Model(CPLEX.Optimizer)

    @variable(sm, 0 <= d1[1:n, 1:n] <= 1)
    @variable(sm, 0 <= d2[1:n, 1:n] <= 2)

    @constraint(sm, sum(d1[i,j] for i in 1:n, j in 1:n) <= T)
    @constraint(sm, sum(d2[i,j] for i in 1:n, j in 1:n) <= T^2)

    @objective(sm, Max, sum((t[i,j] + d1[i,j]*(th[i]+th[j]) + d2[i,j]*th[i]*th[j])*xx[i,j] for i in 1:n, j in 1:n))

    #set_optimizer_attribute(sm, "OutputFlag", 0)

    optimize!(sm)

    if primal_status(sm) == MOI.FEASIBLE_POINT
        return objective_value(sm), value.(d1), value.(d2)
    else
        error("SP infaisable")
    end
end

function greedy_SP(xx)
    t1 = vec([((th[i]+th[j])*xx[i,j],i,j) for i in 1:n, j in 1:n])
    t1 = sort(t1, by=x -> -x[1])
    d1 = zeros(n,n)
    budget1 = T
    i = 1
    while budget1 > 0 && i <= n*n
        x = xx[t1[i][2],t1[i][3]]
        aff = min(1*x,budget1*x)
        d1[t1[i][2],t1[i][3]] = aff
        budget1 -= aff
        i += 1
    end

    t2 = vec([((th[i]*th[j])*xx[i,j],i,j) for i in 1:n, j in 1:n])
    sort(t2, by=x -> -x[1])
    d2 = zeros(n,n)
    budget2 = T*T
    i = 1
    while budget2 > 0 && i <= n*n
        x = xx[t2[i][2],t2[i][3]]
        aff = min(2*x,budget2*x)
        d2[t2[i][2],t2[i][3]] = aff
        budget2 -= aff
        i += 1
    end

    return sum((t[i,j] + d1[i,j]*(th[i]+th[j]) + d2[i,j]*th[i]*th[j])*xx[i,j] for i in 1:n, j in 1:n), d1, d2
end

function branch_and_cut()
    m = Model(CPLEX.Optimizer)

    @variable(m, z, Int)
    @variable(m, x[1:n, 1:n], Bin)
    @variable(m, u[2:n], Int)
    
    @constraint(m, z >= sum(t[i,j]*x[i,j] for i in 1:n, j in 1:n))
    @constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[1,j] for j in 2:n) == sum(x[j,1] for j in 2:n))
    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n, i != j], u[j] - u[i] >= d[i] - C*(1 - x[i,j]))
    @constraint(m, [j in 2:n], u[j] <= C*(1 - x[1,j]))
    
    @objective(m, Min, z)
    
    #set_optimizer_attribute(m, "OutputFlag", 0)
        
    i = 1
    while true
        println("Iteration : ", i)
        optimize!(m)

        if primal_status(m) != MOI.FEASIBLE_POINT
            println("Aïe : ", primal_status(m))
            return
        end
    
        zz = objective_value(m)
        xx = value.(x)
        sz, dd1, dd2 = greedy_SP(xx)

        if sz <= zz
            println("Valeur optimale : ", objective_value(m))
            println("Solution x : ", value.(x))        
            break
        end

        @constraint(m, z >= sum((t[i,j] + dd1[i,j]*(th[i]+th[j]) + dd2[i,j]*th[i]*th[j])*x[i,j] for i in 1:n, j in 1:n))

        i += 1
    end
end

println("Instance : ", file)
branch_and_cut()
