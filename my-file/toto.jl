function main(xu)
a::eltype(xu) = 1.
return a
end

b = main([4.,7])

println("typeof(b) = ", typeof(b))

