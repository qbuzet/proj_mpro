using JuMP
#using CPLEX
using GLPK

function generateBlock(n::Int, block_number::Int)
    blockSize = round(Int, sqrt(n))
    i_start = blockSize * (div(block_number - 1, blockSize) + 1) - blockSize + 1
    j_start = blockSize * (mod(block_number - 1, blockSize) + 1) - blockSize + 1
    return [(i, j) for i in i_start:i_start+blockSize-1, j in j_start:j_start+blockSize-1]
end

function robuste_dualisation()

  