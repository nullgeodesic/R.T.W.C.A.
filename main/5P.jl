"""
Author: Levi Malmström
"""

using Pkg
#Pkg.add("CUDA")
using CUDA
#Pkg.add("StaticArrays")
using StaticArrays
using Base.Threads

include("integrators_5P.jl")
include("tests_5P.jl")


#println(CUDA.versioninfo())

"""
Runs the integration loop for a single pixel.
"""
@inline function integrate_ray!(ray,starting_timestep,colors,colors_freq,::Val{raylength},abs_tol,rel_tol,max_dt_scale,
                    max_steps) where raylength
    #integrate ray
    dt = starting_timestep
    
    #various 'buffers' to drastically reduce memory allocations
    source_vel = @MVector zeros(Float32,4)
    g = @MMatrix zeros(Float32,4,4)
    shared_slope = @MVector zeros(Float32,raylength)
    last_slope = @MVector zeros(Float32,raylength)
    next_slope = @MVector zeros(Float32,raylength)

    #calculate initial derivative
    calc_ray_derivative!(ray,raylength,colors_freq,shared_slope,source_vel,g)
    @inbounds for i in 1:raylength
        last_slope[i] = shared_slope[i]
        next_slope[i] = shared_slope[i]
    end

    #more buffers
    buffer = @MVector zeros(Float32,raylength)
    y = @MVector zeros(Float32,raylength)

    k2 = @MVector zeros(Float32,raylength)
    k3 = @MVector zeros(Float32,raylength)
    k4 = @MVector zeros(Float32,raylength)
    k5 = @MVector zeros(Float32,raylength)
    k6 = @MVector zeros(Float32,raylength)
    k7 = @MVector zeros(Float32,raylength)

    rejected = false
    raytrace=true
    step_count = 0

    
    while raytrace
        dt,rejected = RKDP_Step_w_buffer!(ray,y,last_slope,next_slope,raylength,dt,colors_freq,
                                                         abs_tol, rel_tol,rejected,max_dt_scale,buffer,
                                                         k2,k3,k4,k5,k6,k7,source_vel,g)
        step_count+=1

        
        #my current termination condition
        if calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,
                          max_dt_scale, max_steps,step_count)
            raytrace=false
        end
    end
end


function ray_kernel!(long_ray_matrix,colors,colors_freq,
                     ::Val{raylength},abs_tol,rel_tol,max_dt_scale,max_steps,num_pix) where raylength
    
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x    
    stride = gridDim().x * blockDim().x
    
    for i = index:stride:num_pix
        ray = @MVector zeros(Float32,raylength)
        for j in 1:raylength
            @inbounds ray[j] = long_ray_matrix[i,j]
        end
        
        starting_timestep = -pad_max_dt(ray,max_dt_scale)
        
        integrate_ray!(ray,starting_timestep,colors,colors_freq,Val(raylength),abs_tol,rel_tol,max_dt_scale,
                            max_steps)
        
        
        for j in 1:raylength
            @inbounds long_ray_matrix[i,j] = ray[j]
        end
        
    end
    
    return nothing
end


function color_kernel!(ray_bundle,img_bundle,colors,colors_freq,raylength::Integer,beamsize::Real)
    for i in 1:size(ray_bundle,1)
        @views ray = ray_bundle[i,:]
        skybox_handling!(ray,raylength,colors,colors_freq,beamsize)
        #calculate pixel value
        img_bundle[i] = calc_xyY(ray,colors,colors_freq)
    end
    return img_bundle
end

                       
"""
    gen_image(;camera_pos=[0,0,0,0],camera_dir=[0.0,0.0],speed=0.0::Real,
                        camera_point=[0.0,0.0,0.0],x_pix=40,max_steps=1e3,colors=[400,550,700],
                             returnrays=false,tolerance=1e-4::Real,max_dt_scale=1e-2::Real,fov_hor=85::Real,
                   fov_vert=60::Real)
Generates an image. Multithreaded.
# Arguments
- 'camera_pos': the 4-position of the camera
- 'camera_dir': two angles determining the direction of the camera's 3-velocity relative to the tetrad/FIDO.
- 'speed': the magnitude of the camera's 3-velocity relative to the tetrad/FIDO in units of c.
- 'camera_point': the orientation of the camera relative to it's velocity vector.
- 'x_pix': the size of the bottom of the image, in pixels.
- 'max_steps': the maximum number of timesteps to integrate.
- 'colors': the wavelengths of frequencies to use in the spectrum, in nm.
- 'returnrays': return the rays which were traced.
- 'tolerance': number I'm currently using to set the error tolerance of the integrator.
- 'max_dt_scale': a positive number characterising the maximum size of timesteps. Decrease to shrink the timestep size.
- 'fov_hor': horizontal field of view, in degrees.
- 'fov_vert': vertical field of view, in degrees.
"""
function gen_image(;camera_pos=[0,0,0,0],camera_dir=[0.0,0.0],speed=0.0::Real,
                        camera_point=[0.0,0.0,0.0],x_pix=40,max_steps=1e3,colors=[400,550,700],
                             returnrays=false,tolerance=1e-4::Real,max_dt_scale=1e-2::Real,fov_hor=85::Real,
                   fov_vert=60::Real,print_num_pix=false::Bool)
    #check that camera is in a valid location
    if is_singularity(camera_pos)
        println("Invalid Camera; returning blank image")
        return zeros(RGB{N0f16},2)
    end

    #initialize rays
    ray_matrix, beamsize = initialize_camera(camera_pos,direction=camera_dir,β=speed,pointing=camera_point,
                                   horizontal_pixels=x_pix,colors=colors,fov_hor=fov_hor, fov_vert=fov_vert)
    y_pix=size(ray_matrix,2)
    num_pix=x_pix*y_pix
    if print_num_pix
        println("Number of Pixels ",num_pix)
    end

    n_colors=length(colors)
    raylength = 8+2*n_colors
    long_ray_matrix = reshape(ray_matrix,(num_pix,raylength))
    ray_matrix_shape = size(ray_matrix)
    
    #initialize xyY pixels
    xyY_img = zeros(xyY{Float64},x_pix,ray_matrix_shape[2])
    long_xyY_img = reshape(xyY_img,num_pix)

    f(x)=1/x

    #initialize colors_freq
    colors_freq = f.(colors)

    #make tolerance arrays
    abs_tol=fill(tolerance,raylength)
    rel_tol=fill(tolerance,raylength)



    #transfer/transfrom data to Device
    Cu_ray_matrix = CuArray{Float32}(long_ray_matrix)
    Cu_colors = CuArray{Float32}(colors)
    Cu_colors_freq = CuArray{Float32}(colors_freq)
    Cu_abs_tol = CuArray{Float32}(abs_tol)
    Cu_rel_tol = CuArray{Float32}(rel_tol)

    #code for debugging kernel
    """
    @device_code_llvm raw=true dump_module=true @cuda ray_kernel!(Cu_ray_matrix,Cu_colors,Cu_colors_freq,
                                                   Val(raylength),Cu_abs_tol,Cu_rel_tol,Float32(max_dt_scale),
                                                   Int32(max_steps),Int32(num_pix))
    """
    
    #compile the kernel and make the configuration
    Cu_ray_kernel! = @cuda launch=false ray_kernel!(Cu_ray_matrix,Cu_colors,Cu_colors_freq,
                                                   Val(raylength),Cu_abs_tol,Cu_rel_tol,Float32(max_dt_scale),
                                                   Int32(max_steps),Int32(num_pix))
    config = launch_configuration(Cu_ray_kernel!.fun)
    threads = min(num_pix,config.threads)
    blocks = cld(num_pix,threads)
    println("Threads: ", threads)
    println("Blocks: ", blocks)

    
    #Integrate the rays on the Device
    
    CUDA.@sync begin
        Cu_ray_kernel!(Cu_ray_matrix,Cu_colors,Cu_colors_freq,
                                                   Val(raylength),Cu_abs_tol,Cu_rel_tol,Float32(max_dt_scale),
                                                   Int32(max_steps),Int32(num_pix); threads, blocks)
    end
    
    #Bring back the integrated rays to the host
    long_ray_matrix = Array(Cu_ray_matrix)

    #Calculate the colors of pixels
    chunks = Iterators.partition(1:num_pix,cld(num_pix,100*nthreads()))
    tasks = []
    tasks = map(chunks) do chunk
        @spawn color_kernel!(deepcopy(long_ray_matrix[chunk,:]),deepcopy(long_xyY_img[chunk]),
                            deepcopy(colors),deepcopy(colors_freq),raylength,beamsize)
    end
    finished_bundles = fetch.(tasks)
    index_counter = 1
    for i in 1:length(finished_bundles)
        for j in 1:size(finished_bundles[i],1)
            long_xyY_img[index_counter] = finished_bundles[i][j]
            index_counter += 1
        end
    end

    #transform back into a grid
    ray_matrix=reshape(long_ray_matrix,(x_pix,y_pix,raylength))
    xyY_img=reshape(long_xyY_img,(x_pix,y_pix))

    #check that any rays even hit anything, and return a blank image if they didn't
    max_pixel_val = maxfinite(xyY_img)
    if comp3(max_pixel_val) <= 0
        println("Image blank")
        if returnrays
            return zeros(RGB{N0f16},size(ray_matrix,2),size(ray_matrix,1)),ray_matrix
        else
            return zeros(RGB{N0f16},size(ray_matrix,2),size(ray_matrix,1))
        end
    end
    
    #scale the image
    scaler = scaleminmax(0,comp3(max_pixel_val))
    for i in 1:size(xyY_img,1)
        for j in 1:size(xyY_img,2)
            xyY_img[i,j] = xyY{Float64}(comp1(xyY_img[i,j]),comp2(xyY_img[i,j]),scaler(comp3(xyY_img[i,j])))
        end
    end
    
    
    #convert to rgb and return
    RGB_img = convert.(RGB{N0f16},xyY_img)
    if returnrays
        return transpose(reverse(RGB_img,dims=2)),ray_matrix
    else
        return transpose(reverse(RGB_img,dims=2))
    end
end
