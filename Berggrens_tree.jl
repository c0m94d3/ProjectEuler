function primitive_triples(depth::Int64, max_sum::Int64, list::Bool = false)
    node_0 = Vector{Int64}([3; 4; 5])
    nodes = Dict{Int64, Vector{Vector{Int64}}}([0 => [node_0]])
    list && (ptriples = Vector{Vector{Int64}}([node_0]))
    Barning_matrix_1 = Matrix{Int8}([[1 -2 2]; [2 -1 2]; [2 -2 3]])
    Barning_matrix_2 = Matrix{Int8}([[1 2 2]; [2 1 2]; [2 2 3]])
    Barning_matrix_3 = Matrix{Int8}([[-1 2 2]; [-2 1 2]; [-2 2 3]])
    nodesofar = 0
    count = 0
    while nodesofar < depth
        triples = Vector{Vector{Int64}}()
        for i in nodes[nodesofar]
            vec1 = Barning_matrix_1*i
            vec2 = Barning_matrix_2*i
            vec3 = Barning_matrix_3*i
            push!(triples, vec1, vec2, vec3)
            list && push!(ptriples, vec1, vec2, vec3)
            if max_sum > 0
                if sum(vec1) ≥ max_sum || sum(vec2) ≥ max_sum || sum(vec2) ≥ max_sum
                    nodes[nodesofar + 1] = triples
                    list ? (return ptriples) : (return nodes, count + 1)
                end
            end
            count += 3
        end
        nodes[nodesofar + 1] = triples
        nodesofar += 1
    end
    list ? (return ptriples) : (return nodes, count + 1)
end
function unique_triangle(L::Int64)
    ptriples = primitive_triples(13, L, true)
    count = 0
    for i in ptriples
        if L % sum(i) == 0
            count += 1
        end
        count ≥ 2 && (return false)
    end
    count == 1 ? (return true) : (return nothing)
end
function unique_tr(N::Int64)
    counter = 0
    for l in 2:2:N
        if unique_triangle(l) == true
            counter += 1
        end
    end
    return (counter)
end
@time println(primitive_triples(2, 0, false))