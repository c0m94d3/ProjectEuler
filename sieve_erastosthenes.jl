#using BenchmarkTools
function sieve(N::Int64)
    arr = repeat([true], N)
    for i ∈ 2:N
        if arr[i]
            for j ∈ i+i:i:N
                arr[j] = false
            end
        end
    end
    return [i for i ∈ 2:N if arr[i]]
end

#time = @elapsed sieve(100000000)
#println("Time Elapsed: $(time)s")