"""
v0.3.4
November 10 2025
Author: Levi Malmström
"""

using LinearAlgebra
using Pkg
#Pkg.add("Plots")
using Plots
#Pkg.add("PlotlyBase")
#Pkg.add("PlotlyKaleido")
plotly()

#Pkg.add("Interpolations")
using Interpolations
#Pkg.add("Integrals")
using Integrals
#Pkg.add("Images")
using Images
#Pkg.add("ImageMagick")
using ImageMagick
#Pkg.add("ImageView")
using ImageView
#Pkg.add("FileIO")
using FileIO
#Pkg.add("ProfileView")
#using ProfileView
#Pkg.add("ProfileCanvas")
using ProfileCanvas
#Pkg.add("StaticArrays")
#using StaticArrays
#Pkg.add("LoopVectorization")
using LoopVectorization
using Profile
using Base.Threads

#CONSTANTS
#Planck Constant in J/Hz
const h = 6.62607015e-34
#Reduced Planck Constant in J*s
const hbar = h/(2*pi)
#Speed of light in m/s
const c = 299792458
#Boltzman Constant in J/K
const k_B = 1.380649e-23
#How many meters corespond to one unit of the map
const map_scale = 1
#Protects from dividing by zero in certain situations
const no_div_zero = 1e-24
#Color match function fit values
const cie_matrix=[
[0.362 1.056 -0.065 0.821 0.286 1.217 0.681];
[442.0 599.8 501.1 568.8 530.9 437.0 459.0];
[0.0624 0.0264 0.0490 0.0213 0.0613 0.0845 0.0385];
[0.0374 0.0323 0.0382 0.0247 0.0322 0.0278 0.0725]]
#Dormand-Prince Butcher tableau
const DP_Butcher= [
[0 0 0 0 0 0 0 0];
[1/5 1/5 0 0 0 0 0 0];
[3/10 3/40 9/40 0 0 0 0 0];
[4/5 44/45 -56/15 32/9 0 0 0 0];
[8/9 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0];
[1 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0];
[1 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]]

include("gen_geometry.jl")
include("color_theory.jl")
include("integrators.jl")
include("Spherical_Minkowski_Sphere.source/source_main.jl")
include("tests.jl")


function Heaviside(x)
    if x >= 0
        return 1
    else
        return 0
    end
end


"""
    initialize_camera(position;direction=[0.0,0.0],beta=0.0::Real,pointing=[0.0,0.0,0.0],
       xpixels=1024::Integer,fov_hor=85::Real,fov_vert=60::Real,colors=[400,550,700])
Initialize the camera/rays from user input.
# Arguments
- 'position': the 4-position of the camera
- 'direction': two angles determining the direction of the camera's 3-velocity relative to the tetrad/FIDO.
- 'beta': the magnitude of the camera's 3-velocity relative to the tetrad/FIDO in units of c.
- 'pointing': the orientation of the camera relative to it's velocity vector.
- 'horizontal_pixels=1024': the size of the bottom of the image, in pixels.
- 'fov_hor=85': the width of the image/field of view, in degrees.
- 'fov_vert=60': the target height of the image/field of view, in degrees.
- 'colors=[400,550,700]': the wavelengths in nm that will be used to create the image.
"""
function initialize_camera(position::Vector; direction=[0.0,0.0], beta=0.0::Real, pointing=[0.0,0.0,0.0],
                           horizontal_pixels=1024::Integer, fov_hor=85::Real, fov_vert=60::Real, colors=[400,550,700])
    #convert to radians
    fov_hor=π*fov_hor/180
    fov_vert=π*fov_vert/180
    #calculate canvas size
    W=tan(fov_hor/2)
    H=tan(fov_vert/2)
    #calculate number of vertical pixels
    vertical_pixels=round(Int,horizontal_pixels*H/W)
    #calculate pixel density
    rho_hor=2*W/(horizontal_pixels+1)
    rho_vert=2*H/(vertical_pixels+1)
    
    #initialize blank rays as state vectors (4 position values, 4 velocity values, 2 optical values for each color)
    ray_state_length=8+2*length(colors)
    S=zeros(horizontal_pixels,vertical_pixels,ray_state_length)
    #now actually fill in proper values into those spots!
    for i = 1:horizontal_pixels
        w=rho_hor*(i-0.5)-W
        for j = 1:vertical_pixels
            for k =1:4
                S[i,j,k]=position[k]
            end
            
            h=rho_vert*(j-0.5)-H
            θ=atan(sqrt(w^2 + h^2))

            #phi=sign(h)*acos(w/sqrt(w^2 + h^2))
            if w>0
                ϕ =atan(h/w)
            elseif w<0 && h>=0
                ϕ=atan(h/w)+ π
            elseif w<0 && h<0
                ϕ=atan(h/w)-π
            elseif w==0 && h>0
                ϕ=π/2
            elseif w==0 && h<0
                ϕ=-π/2
            else
                ϕ=0
            end
      
            S[i,j,5]=1
            S[i,j,6]=-sin(θ)*cos(ϕ)
            S[i,j,7]=-sin(θ)*sin(ϕ)
            S[i,j,8]=-cos(θ)
            #geometrically, we have projected the pixels from the canvas at z=1 onto a sphere
            #but the points actually are velocity vectors, so we need to negate them,
            #since the rays are traveling into the camera, not out.
        end
    end

    #Matrix to rotate the camera
    R=gen_intrinsic_rotation_matrix(pointing)
    #Matrix to boost to the FIDO frame
    L=lorentz_boost_z(-beta)
    #Matrix to rotate coordinate axes to align with FIDO frame
    R2=gen_final_rotation_matrix(direction)
    #Matrix to convert to global coordinates
    e=calc_inv_vierbein(position)
    for i = 1:horizontal_pixels
        for j = 1:vertical_pixels
            S[i,j,5:8]=e*R2*L*R*S[i,j,5:8]
        end
    end
    return S
end


function ray_kernel(ray,starting_timestep,tolerance,colors,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale,
                    max_steps)
    #integrate ray
    dt=starting_timestep
    raytrace=true
    shared_slope=calc_ray_derivative(ray,raylength,colors_freq)
    buffer=similar(shared_slope)
    deriv_buffer=similar(shared_slope)
    rejected=false
    step_count=0
    
    while raytrace
        ray,shared_slope,dt,rejected,buffer=RKDP_Step_w_buffer(ray,shared_slope,raylength,dt,colors_freq,
                                                     abs_tol, rel_tol,rejected,max_dt_scale,buffer)
        step_count+=1

        #catch NaN's
        if find_NaN(ray)
            println("NaN Warning, stopping ray ",threadid(), " ", ray[1:8])
            raytrace=false
        end
        
        #my current termination condition
        @views if step_count >= max_steps || minimum(ray[9:2:end]) >= -log(tolerance) || ray[2] > 10
            if step_count >= max_steps
                println("Max Steps",ray[1:8])
            end

            raytrace=false
        end
    end
    
    #calculate pixel value
    xyY_pix=calc_xyY(ray,colors,colors_freq)
    return ray,xyY_pix
end

function iterate_on_pix_range(ray_bundle,img_bundle,tolerance::Real,colors,colors_freq,
                              raylength::Integer,abs_tol,rel_tol,max_dt_scale::Real,max_steps::Real,
                              start_time::UInt64)
    for i in 1:size(ray_bundle,1)
        ray_bundle[i,:],img_bundle[i]=ray_kernel(ray_bundle[i,:],
                                                     -pad_max_dt(ray_bundle[i,1:8],max_dt_scale),
                                                     tolerance,colors,colors_freq,raylength,abs_tol,rel_tol,
                                                     max_dt_scale,max_steps)
    end
    #println("Elapsed time (ns): ",time_ns()-start_time, "; Thread ID: ",threadid())
    #if find_NaN(img_bundle)
        #println(find_NaN(ray_bundle)," ",threadid()," ",ray_bundle[1,1:8])
    #end
    
    return ray_bundle,img_bundle
end


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
    ray_matrix=initialize_camera(camera_pos,direction=camera_dir,beta=speed,pointing=camera_point,
                                 horizontal_pixels=x_pix,colors=colors,fov_hor=fov_hor, fov_vert=fov_vert)
    #initialize xyY pixels
    #xyY_img=Array{xyY{Float64}}(undef,size(ray_matrix,1),size(ray_matrix,2))
    xyY_img=zeros(xyY{Float64},x_pix,size(ray_matrix,2))
    
    n_colors=length(colors)
    raylength=8+2*n_colors
    f(x)=c*1e9/x

    #initialize colors_freq
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
    chunks = Iterators.partition(1:num_pix,cld(num_pix,nthreads()))
    tasks = []
    tasks = map(chunks) do chunk
        @spawn iterate_on_pix_range(deepcopy(long_ray_matrix[chunk,:]),deepcopy(long_xyY_img[chunk]),
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
