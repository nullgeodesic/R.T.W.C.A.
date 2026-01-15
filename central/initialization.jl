"""
v0.4.1
January 15 2026
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
using Profile
#Pkg.add("Cthulhu")
#using Cthulhu


#load general modules
include("gen_geometry.jl")
include("color_theory.jl")
include("special_functions.jl")
include("tests.jl")
include("64BitConstants.jl")
include("32BitConstants.jl")


if mode == "LttM"
    #CPU integration setup
    include("LttM.jl")
    #π in Float64
    const π = Float64(π)
    #Float64 constants
    using .BitConsts64
elseif mode == "5P"
    #GPU integration setup
    include("5P.jl")
    #π in Float32
    const π = Float32(π)
    #Float32 constants
    using .BitConsts32
end


"""
Convenience functions for setup.
"""
function initialize_world(input_string::String)
    input_string = lowercase(input_string)
    if input_string == "schwarzschild_t,r,θ,ϕ_simpledisk"
        include("Schwarzschild_Simple_Disk_Spherical.source/source_main.jl")
        
        position = [0,25,π/3,0]
        
        direction = [π/2,0]
        
        pointing = [0.0,0.0,0.0]
        
        β = 0.0
    elseif input_string == "minkowski_t,r,θ,ϕ_sphere"
        include("Minkowski_Sphere_Spherical.source/source_main.jl")
        
        position = [0,4,π/2,0]
        
        direction = [π/2,0]
        
        pointing = [0.0,0.0,0.0]
        
        β = 0.0
    elseif input_string == "minkowski_t,x,y,z_sphere"
        include("Minkowski_Sphere_Cartesian.source/source_main.jl")

        position = [0.0,0.0,0.0,-4]

        direction = [0.0,0.0]

        pointing = [0.0,0.0,0.0]

        β = 0.0
    else
        println("Type source startup file name e.g. Schwarzschild_Simple_Disk_Spherical.source/source_main.jl")
        source_string = readline()
        include(source_string)
        
        println("Type source starting 4-position e.g. [0,25,π/3,0]")
        position = Float64.(eval.(Meta.parse.(split(strip(readline(),['[',']',' ']), ","))))
        
        println("Type direction-angles of camera's motion relative to the local tetrad at it's location
 e.g. [π/2,0]")
        direction = Float64.(eval.(Meta.parse.(split(strip(readline(),['[',']',' ']), ","))))
        
        println("Type rotation-angles for camera's orientation relative to it's direction of motion
            e.g. [0,0,0]")
        pointing = Float64.(eval.(Meta.parse.(split(strip(readline(),['[',']',' ']), ","))))
        
        println("Type the speed of the camera in units of c; -1 < β < 1; negative values mean the camera is actually
            moving opposite to the direction indicated by the direction-angles")
        β = Float64(eval(Meta.parse(strip(readline(),['[',']',' ']))))
    end

    return position,direction,pointing,β
end

        
"""
    initialize_camera(position;direction=[0.0,0.0],β=0.0::Real,pointing=[0.0,0.0,0.0],
       xpixels=1024::Integer,fov_hor=85::Real,fov_vert=60::Real,colors=[400,550,700])
Initialize the camera/rays from user input.
# Arguments
- 'position': the 4-position of the camera
- 'direction': two angles determining the direction of the camera's 3-velocity relative to the tetrad/FIDO.
- 'β': the magnitude of the camera's 3-velocity relative to the tetrad/FIDO in units of c.
- 'pointing': the orientation of the camera relative to it's velocity vector.
- 'horizontal_pixels=1024': the size of the bottom of the image, in pixels.
- 'fov_hor=85': the width of the image/field of view, in degrees.
- 'fov_vert=60': the target height of the image/field of view, in degrees.
- 'colors=[400,550,700]': the wavelengths in nm that will be used to create the image.
"""
function initialize_camera(position::Vector; direction=[0.0,0.0], β=0.0::Real, pointing=[0.0,0.0,0.0],
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
    L=lorentz_boost_z(-β)
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
