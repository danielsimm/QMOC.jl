function async_task()
    println("spawned task")
    sleep(1)
    println("task done")
    return nothing
end

function wrapper()
    println("main function called")
    Threads.@spawn async_task()
    sleep(1)
    println("main function done")
    return nothing
end


function wrapper_waiting()
    println("main function called")
    task = Threads.@spawn async_task()
    sleep(1)
    println("main function done")
    return nothing
    wait(task)
end

wrapper()