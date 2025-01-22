using JuMP
using Gurobi

file = "data/n_6-euclidean_true"
include(file)

function solve_SP(xx)
    sm = Model(Gurobi.Optimizer)

    @variable(sm, 0 <= d1[1:n, 1:n] <= 1)
    @variable(sm, 0 <= d2[1:n, 1:n] <= 2)

    @constraint(sm, sum(d1[i,j] for i in 1:n, j in 1:n) <= T)
    @constraint(sm, sum(d2[i,j] for i in 1:n, j in 1:n) <= T^2)

    @objective(sm, Max, sum((t[i,j] + d1[i,j]*(th[i]+th[j]) + d2[i,j]*th[i]*th[j])*xx[i,j] for i in 1:n, j in 1:n))

    set_optimizer_attribute(sm, "OutputFlag", 0)

    optimize!(sm)

    if primal_status(sm) == MOI.FEASIBLE_POINT
        return objective_value(sm), value.(d1), value.(d2)
    else
        error("SP infaisable")
    end
end

counter = 0

function my_callback_function(cb_data)
    global counter
    counter += 1

    status = callback_node_status(cb_data, m)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        return
    elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER

    else
        @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN

    end

    zz = callback_value(cb_data, z)
    xx = [callback_value(cb_data, x[i, j]) for i in 1:n, j in 1:n]

    println("Callback n°", counter, ", z* : ", zz)

    sz, dd1, dd2 = solve_SP(xx)
    
    if sz > zz
        ncon = @build_constraint(z >= sum((t[i,j] + dd1[i,j]*(th[i]+th[j]) + dd2[i,j]*th[i]*th[j])*x[i,j] for i in 1:n, j in 1:n))
        MOI.submit(m, MOI.LazyConstraint(cb_data), ncon)
    end
    return
end

m = Model(Gurobi.Optimizer)

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

set_optimizer_attribute(m, "OutputFlag", 0)
set_attribute(m, MOI.LazyConstraintCallback(), my_callback_function)

optimize!(m)

if primal_status(m) == MOI.FEASIBLE_POINT
    println("Valeur optimale : ", objective_value(m))
    println("Solution x : ", value.(x))
else
    println("Aïe : ", primal_status(m))
end