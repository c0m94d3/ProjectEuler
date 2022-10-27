#module partitions
using BenchmarkTools
using DataStructures
using DSP
#export ZS1, ZS2, partitions_recurrence, partitions_generatingfn
#end
#=
ZS1 and ZS2 based on ALGORITHMS explained in :
FAST ALGORITHMS FOR GENERATING
INTEGER PARTITIONS
ANTOINE ZOGHBIU and IVAN STOJMENOVIC*
Bell Northern Research, P.O . Box 35 I I , Station C , Mail Stop 091 ,
Ottawa, Ontario KIY 4H7, Canada;
b Computer Science Department, SITE, University of Ottawa,
Ottawa. Ontario KIN 6N5. Canada
=#

function ZS1(n::Int64)
    partitions = Array{Array{Int64, 1}, 1}()
    x_i = Dict{Int64, Int64}([i=>1 for i in 1:n])
    x_i[1] = n
    m = 1
    h = 1
    push!(partitions, Vector{Int64}([x_i[1]]))
    while x_i[1] ≠ 1
        if x_i[h] == 2
            m += 1
            x_i[h] = 1
            h = h - 1
        else
            r = x_i[h] - 1
            t = m - h + 1
            x_i[h] = r
            while t ≥ r
                h = h + 1
                x_i[h] = r
                t = t - r
            end
            if t == 0
                m = h
            else
                m = h + 1
                if t > 1
                    h = h + 1
                    x_i[h] = t 
                end
            end
        end
        push!(partitions, Array{Int64, 1}([x_i[j] for j∈1:m]))
    end
    return partitions
end
function ZS2(n::Int64)
    partitions = Array{Array{Int64, 1}, 1}()
    x_i = Dict{Int64, Int64}([i=>1 for i in 1:n])
    push!(partitions, [1 for j in 1:n])
    x_i[0] = -1
    x_i[1] = 2
    h = 1
    m = n - 1
    push!(partitions, [x_i[j] for j in 1:m])
    while x_i[1] ≠ n
        if m - h > 1
            h = h + 1
            x_i[h] = 2
            m = m - 1
        else
            j = m - 2
            while x_i[j] == x_i[m-1]
                x_i[j] = 1
                j = j - 1
            end
            h = j + 1
            x_i[h] = x_i[m-1] + 1
            r = x_i[m] + x_i[m-1]*(m-h-1)
            x_i[m] = 1
            if m - h > 1
                x_i[m-1] = 1
            end
            m = h + r - 1
        end
        push!(partitions, Array{Int64, 1}([x_i[j] for j∈1:m]))
    end
    return partitions
end
function partitions_recurrence(N::Int64, k::Int64, dict::Bool=false, ordered::Bool=false)
#= Generates p(n) (unrestricted partition function) ∀ n ∈ 1 to N, returns p(k)
Based on the recurrence relation: 
p(n) = p(n-1) + p(n-2) + p(n-5) + p(n-7) + ...
which 
=#
    generalized_pentnums = k -> Tuple{Int64, Int64}([((3k^2-k)//2).num, ((3k^2+k)//2).num])
    if ordered
        cached_p = OrderedDict{Int64, BigInt}([0 => 1, 1 => 1, 2 => 2])
    else
        cached_p = Dict{Int64, BigInt}([0 => 1, 1 => 1, 2 => 2])
    end
    for j in 3:N
        p_j::BigInt = 0
        i = 1
        Λ = Tuple{Function, Function}((+, -))
        while true
            exponents = generalized_pentnums(i)
            if exponents[1] > j
#                println("Breaking, $j-$(exponents[1]) < 0")
                break
            else
#                println("p($j) $(Λ[1])= p($(j-exponents[1]))")
                p_j = Λ[1](p_j, cached_p[j-exponents[1]])
            end
            if exponents[2] > j
#                println("Breaking, $j-$(exponents[2]) < 0")
                break
            else
#                println("p($j) $(Λ[1])= p($(j-exponents[2]))")
                p_j = Λ[1](p_j, cached_p[j-exponents[2]])
            end
            Λ = reverse(Λ)
            i += 1
        end
        cached_p[j] = p_j
    end
    dict && return cached_p;
    return cached_p[k];
end
function fn_series(n::Int64, length::Int64) :: Vector{Int8}
#=
Converts a generating function of the form 1/(1-x^n) to the formal power series it represents.
1/(1-x^i) = 1 + x^i + x^2i + x^3i + ... (length)
=#
    series = zeros(Int8, length + 1)
    series[1] = Int8(1)
    i = 1
    while i < (length + 1)//n
        series[i*n + 1] = Int8(1)
        i += 1
    end
    return series
end
function partitions_generatingfn(N::Int64, denominations::Array{Int64, 1})
    if length(denominations) > 0
        return reduce((i, j) -> conv(i, j), Vector{Vector{Int8}}([fn_series(j, N) for j in denominations]))[N + 1]
    else
        return reduce((i, j) -> conv(i, j), Vector{Vector{Int8}}([fn_series(j, N) for j in 1:N]))[N + 1]
    end
end
#@time println(partitions_generatingfn(10000, [1, 2, 5, 10, 20, 50, 100, 200]))
#@time println(partitions_recurrence(1000, 1000))