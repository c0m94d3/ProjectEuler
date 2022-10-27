using ResumableFunctions
using BenchmarkTools
function partial_quotients(N::Int64)
    issquare(n) = isqrt(n)^2 == n
    if issquare(N)
        return Array([N, Array{Int8, 1}()])
    end
    quotients = Array{Union{Int64, Array{Int64, 1}}, 1}([isqrt(N), Array{Int64, 1}()])
    m = 0
    d = 1
    a = isqrt(N)
    count = 0
    while a ≠ 2isqrt(N)
        count += 1
        mn = d*a-m
        dn = (N - mn^2)/d
        an = floor((isqrt(N)+mn)/dn)
        append!(quotients[2], an)
        m, d, a = mn, dn, an
    end
    return quotients
end
@resumable function convergents(N::Int64)# Q::Int32
    partials = Iterators.cycle(partial_quotients(N)[2])
#    convergent = Array{Tuple{Int64, Int64}, 1}()
    nume_0::Int128 = 1
    denom_0::Int128 = 0
    nume_::Int128 = isqrt(N)
    denom_::Int128 = 1
#    push!(convergent, (nume_, denom_))
#    count = 1
    for i in partials
        while true #count < Q + 1
            nume1, denom1 = nume_, denom_
            nume_, denom_ = nume_*i+nume_0, denom_*i+denom_0
            nume_0, denom_0 = nume1, denom1
#            push!(convergent, (nume_, denom_))
            @yield Tuple{Int128, Int128}((nume_, denom_))
#            count += 1
            break
        end
        #=
        if count==Q
            break
        end
        =#
    end
#    return convergent
end
function pell_solutions(D::Int64) #minimal solution for x^2-Dy^2=1; (x, y) ∈ Convergents(√D)
    for i in convergents(D)
        if i[1]^2-D*i[2]^2 == 1
            return i
        end
    end
end
function squarefree_ints(start::Int64, stop::Int64)
    squarefree = Array{Int64, 1}()
    issquare(n) = isqrt(n)^2 == n
    for i ∈ start:stop
        if !issquare(i)
            append!(squarefree, i)
        end
    end
    return squarefree
end
function result()
    nonsquares = squarefree_ints(1, 1000)
    sols::Array{Tuple{Int128, Int128}, 1} = [pell_solutions(i) for i in nonsquares]
    max_::Int128 = 0
    for (count, i) in enumerate(sols)
        if i[1]>max_
            max_ = i[1]
            return max_, nonsquares[count]
#            println("x = $max_\t D = $(nonsquares[count])")
        end
    end
end
@time println(result())