using Primes
function multiplicative_order(B::Int64, N::Int64)
    if N == 1
        return 0
    end
    if gcd(B, N) != 1
        return -1
    end
    i = 1
    j = 1
    while i < N
        j = (j * B) % N
        j == 1 && (return i)
        i += 1
    end
end
function period(N::Int64)
    if gcd(N, 10) == 1
        return multiplicative_order(10, N)
    else
        otherfactors = prodfactors(filter(i -> i.first∉[2, 5], factor(N)))
        return multiplicative_order(10, otherfactors)
    end
end
function solution()
    max_p = 0
    max_d = 0
    for i in 1:1000
        p = period(i)
        if p > max_p
            max_p = p
            max_d = i
        end
    end
    return max_p, max_d
end
@time println(solution())