println("Threads: ",Threads.nthreads())

function set_value(c, i)
    offset = (i-1)*3
    c[offset+1:offset+2] = [(i-1)*2+1, (i-1)*2+2]
end
function set_minus(c, i)
    offset = (i-1)*3+2
    c[offset+1] = -1
end

function set_value_view(c_block, arg)
    c_block[1:2] = [(arg-1)*2+1, (arg-1)*2+2]
end
function set_minus_view(c_block, arg)
    c_block[1] = -1
end

function test_loop()
    N = 5
    c = Vector{Float64}(undef, N*3)
    args = 1:N

    Threads.@threads for i=1:N
        #set_value(c, i)
        #set_minus(c, i)
        set_value_view(@view(c[(i-1)*3+1:(i-1)*3+2]), args[i])
        set_minus_view(@view(c[(i-1)*3+2+1:(i-1)*3+2+1]), args[i])
    end

    println(c)
end