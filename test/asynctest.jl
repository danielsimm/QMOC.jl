function async_task_1()
    println("spawned task A")
    sleep(2)
    println("task done")
    return nothing
end

function async_task_2()
    println("spawned task B")
    sleep(3)
    println("task done")
    return nothing
end

function wrapper()
    for i in 1:10
        @sync begin
            println("run $(i)")
            Threads.@spawn async_task_1()
            Threads.@spawn async_task_2()
        end
    end
    
    return nothing
end


wrapper() #this actually works.