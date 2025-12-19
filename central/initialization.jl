"""
v0.3.6
December 18 2025
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

#load general modules
include("gen_geometry.jl")
include("color_theory.jl")
include("integrators.jl")
include("tests.jl")
include("special_functions.jl")
include("control_multithread.jl")

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
