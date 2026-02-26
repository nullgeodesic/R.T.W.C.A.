"""
v0.4.0
Jan 1 2025
Author: Levi Malmström
"""

function test_kernel(gpu_array,::Val{raylength}) where raylength
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    ray = @MVector zeros(Float32,raylength)

    raysize = Size(ray)
    len = raysize[1]
    for i in 1:len
        ray[i] = gpu_array[index,i]
    end
    @cuprintln(index," ",ray[1]," ",ray[2]," ",ray[3]," ",ray[4])
    return nothing
end


function test_test_kernel(npix::Int,raylength::Int)
    A = ones(Float32,(npix,raylength))
    println(A)
    cu_A = cu(A)

    cu_test_kernel = @cuda launch=false test_kernel(cu_A,Val(raylength))
    config = launch_configuration(cu_test_kernel.fun)
    threads = threads = min(size(A,1),config.threads)
    blocks = blocks = cld(size(A,1),threads)

    CUDA.@sync begin
               cu_test_kernel(cu_A; threads, blocks)
    end
end
