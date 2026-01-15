"""
v0.4.1
January 15 2026
Author: Levi Malmström
"""

include("integrators_LttM.jl")

using Base.Threads

"""
Runs the integration loop for a single pixel, then calculates it's xyY colorspace value.
"""
function integrate_ray(ray,starting_timestep,tolerance,colors,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale,
                    max_steps)
    #integrate ray
    dt=starting_timestep

    #various 'buffers' to drastically reduce memory allocations
    source_vel = [0.0,0.0,0.0,0.0]
    g = Matrix{Float64}(I,4,4)
    shared_slope = Array{Float64}(undef,raylength)
    calc_ray_derivative!(ray,raylength,colors_freq,shared_slope,source_vel,g)
    last_slope = copy(shared_slope)
    next_slope = copy(shared_slope)
    buffer = similar(shared_slope)
    y = similar(shared_slope)

    k2 = similar(shared_slope)
    k3 = similar(shared_slope)
    k4 = similar(shared_slope)
    k5 = similar(shared_slope)
    k6 = similar(shared_slope)
    k7 = similar(shared_slope)
    
    raytrace=true
    rejected=false
    step_count=0
    
    while raytrace
        dt,rejected=RKDP_Step_w_buffer!(ray,y,last_slope,next_slope,raylength,dt,colors_freq,
                                                         abs_tol, rel_tol,rejected,max_dt_scale,buffer,
                                                         k2,k3,k4,k5,k6,k7,source_vel,g)
        step_count+=1

        #catch NaN's
        if find_NaN(ray)
            println("NaN Warning, stopping ray ",threadid(), " ", ray[1:8])
            raytrace=false
        end

        #my current termination condition
        if calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,
                          max_dt_scale, max_steps,step_count)
            """
            if step_count >= max_steps
                println("Max Steps",ray[1:8])
            end
            """
            raytrace=false
        end
    end
    
    #calculate pixel value
    xyY_pix=calc_xyY(ray,colors,colors_freq)
    return ray,xyY_pix
end


"""
Solves the rays given to it, feeding them into integrate_ray sequentialy.
"""
function ray_kernel(ray_bundle,img_bundle,tolerance::Real,colors,colors_freq,
                              raylength::Integer,abs_tol,rel_tol,max_dt_scale::Real,max_steps::Real,
                              start_time::UInt64)
    for i in 1:size(ray_bundle,1)
        ray_bundle[i,:],img_bundle[i]=integrate_ray(ray_bundle[i,:],
                                                     -pad_max_dt(ray_bundle[i,1:8],max_dt_scale),
                                                     tolerance,colors,colors_freq,raylength,abs_tol,rel_tol,
                                                     max_dt_scale,max_steps)
    end
    #println("Elapsed time (ns): ",time_ns()-start_time, "; Thread ID: ",threadid())
    
    return ray_bundle,img_bundle
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
                   fov_vert=60::Real)
    #check that camera is in a valid location
    if is_singularity(camera_pos)
        println("Invalid Camera; returning blank image")
        return zeros(RGB{N0f16},2)
    end

    #initialize rays
    ray_matrix=initialize_camera(camera_pos,direction=camera_dir,β=speed,pointing=camera_point,
                                 horizontal_pixels=x_pix,colors=colors,fov_hor=fov_hor, fov_vert=fov_vert)
    #initialize xyY pixels
    #xyY_img=Array{xyY{Float64}}(undef,size(ray_matrix,1),size(ray_matrix,2))
    xyY_img=zeros(xyY{Float64},x_pix,size(ray_matrix,2))
    
    n_colors=length(colors)
    raylength=8+2*n_colors

    #initialize colors_freq (calculate frequency in units of nm^-1)
    f(x)=1/x
    colors_freq=f.(colors)

    #make tolerance arrays
    abs_tol=fill(tolerance,raylength)
    rel_tol=fill(tolerance,raylength)

    y_pix=size(ray_matrix,2)
    num_pix=x_pix*y_pix
    println("num_pix= ",num_pix)

    long_ray_matrix = reshape(ray_matrix,(num_pix,raylength))
    long_xyY_img = reshape(xyY_img,num_pix)
    
    start_time =time_ns()

    #Main Loop
    chunks = Iterators.partition(1:num_pix,cld(num_pix,100*nthreads()))
    tasks = []
    tasks = map(chunks) do chunk
        @spawn ray_kernel(deepcopy(long_ray_matrix[chunk,:]),deepcopy(long_xyY_img[chunk]),
                                    deepcopy(tolerance),deepcopy(colors),deepcopy(colors_freq),
                                    deepcopy(raylength),deepcopy(abs_tol),deepcopy(rel_tol),
                                    deepcopy(max_dt_scale),deepcopy(max_steps),deepcopy(start_time))
    end
    finished_bundles=fetch.(tasks)
    
    index_counter=1
    for i in 1:length(finished_bundles)
        for j in 1:size(finished_bundles[i][1],1)
            long_ray_matrix[index_counter,:]=finished_bundles[i][1][j,:]
            long_xyY_img[index_counter]=finished_bundles[i][2][j]
            index_counter+=1
        end
    end

    ray_matrix=reshape(long_ray_matrix,(x_pix,y_pix,raylength))
    xyY_img=reshape(long_xyY_img,(x_pix,y_pix))
    #println(find_NaN(xyY_img))
    
    #check that any rays even hit anything, and return a blank image if they didn't
    max_pixel_val=maxfinite(xyY_img)
    if comp3(max_pixel_val)<=0
        println("Image blank")
        if returnrays
            return zeros(RGB{N0f16},size(ray_matrix,2),size(ray_matrix,1)),ray_matrix
        else
            return zeros(RGB{N0f16},size(ray_matrix,2),size(ray_matrix,1))
        end
    end
    
    #scale the image
    scaler=scaleminmax(0,comp3(max_pixel_val))
    for i in 1:size(xyY_img,1)
        for j in 1:size(xyY_img,2)
            xyY_img[i,j]=xyY{Float64}(comp1(xyY_img[i,j]),comp2(xyY_img[i,j]),scaler(comp3(xyY_img[i,j])))
        end
    end
    
    
    #convert to rgb and return
    RGB_img=convert.(RGB{N0f16},xyY_img)
    if returnrays
        return transpose(reverse(RGB_img,dims=2)),ray_matrix
    else
        return transpose(reverse(RGB_img,dims=2))
    end
end
