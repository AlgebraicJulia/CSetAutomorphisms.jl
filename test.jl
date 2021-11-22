f() = begin
    sleep(3)
    return 1
end

function timeout(f::Function, n::Int)
    timer = Timer(n)
    while isopen(timer)
        yield()
        return f()
    end
    return nothing
end

println(timeout(f, 2))