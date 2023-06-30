using Revise
function main(xu::Vector{<:Real})
println("xu = ", xu)
a::eltype(xu) = 1
println("typeof(a) = ", typeof(a))
return a
end

c = [4,7]
println("b = main(c)")
println("typeof(c) = ", typeof(c))
b = main(c)
println("typeof(b) = ", typeof(b))

c = [4.,7]
println("b = main(c)")
println("typeof(c) = ", typeof(c))
b = main(c)
println("typeof(b) = ", typeof(b))

#= include("my-file/toto.jl")
b = main(c)
typeof(c) = Vector{Int64}
xu = [4, 7]
typeof(a) = Int64
typeof(b) = Int64
b = main(c)
typeof(c) = Vector{Float64}
xu = [4.0, 7.0]
typeof(a) = Float64
typeof(b) = Float64
=#

