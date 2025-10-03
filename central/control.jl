"""
v0.3.1
October 1 2025
Author: Levi Malmström
"""

using LinearAlgebra
using Pkg
#Pkg.add("Plots")
using Plots
#Pkg.add("Interpolations")
using Interpolations
#Pkg.add("Integrals")
using Integrals
#Pkg.add("Images")
using Images
#Pkg.add("ImageMagick")
using ImageMagick
using ImageView
using FileIO

#CONSTANTS
#Planck Constant in J/Hz
const h=6.62607015e-34
#Reduced Planck Constant in J*s
const hbar=h/(2*pi)
#Speed of light in m/s
const c=299792458
#Boltzman Constant in J/K
const k_B=1.380649e-23
#How many meters corespond to one unit of the map
const map_scale=1
#Color match function fit values
const cie_matrix=[
[0.362 1.056 -0.065 0.821 0.286 1.217 0.681];
[442.0 599.8 501.1 568.8 530.9 437.0 459.0];
[0.0624 0.0264 0.0490 0.0213 0.0613 0.0845 0.0385];
[0.0374 0.0323 0.0382 0.0247 0.0322 0.0278 0.0725]]
#Dormand-Prince Butcher tableau
const DP_Butcher=[
[0 0 0 0 0 0 0 0];
[1/5 1/5 0 0 0 0 0 0];
[3/10 3/40 9/40 0 0 0 0 0];
[4/5 44/45 -56/15 32/9 0 0 0 0];
[8/9 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0];
[1 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0];
[1 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]]

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
function initialize_camera(position; direction=[0.0,0.0], beta=0.0::Real, pointing=[0.0,0.0,0.0],
                           horizontal_pixels=1024::Integer, fov_hor=85::Real, fov_vert=60::Real, colors=[400,550,700])
    #convert to radians
    fov_hor=pi*fov_hor/180
    fov_vert=pi*fov_vert/180
    #calculate canvas size
    W=tan(fov_hor/2)
    H=tan(fov_vert/2)
    #calculate number of vertical pixels
    vertical_pixels=round(Int,horizontal_pixels*H/W)
    #calculate pixel density
    rho_hor=2*W/(horizontal_pixels+1)
    rho_vert=2*H/(vertical_pixels+1)
    
    #initialize blank rays as state vectors (4 position values, 4 velocity values, 2 optical values for each color
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
            theta=atan(sqrt(w^2 + h^2))

            #phi=sign(h)*acos(w/sqrt(w^2 + h^2))
            if w>0
               phi=atan(h/w)
            elseif w<0 && h>=0
                phi=atan(h/w)+pi
            elseif w<0 && h<0
                phi=atan(h/w)-pi
            elseif w==0 && h>0
                phi=pi/2
            elseif w==0 && h<0
                phi=-pi/2
            else
                phi=0
            end
      
            S[i,j,5]=1
            S[i,j,6]=-sin(theta)*cos(phi)
            S[i,j,7]=-sin(theta)*sin(phi)
            S[i,j,8]=-cos(theta)
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
    e=calc_local_vierbein(position)
    for i = 1:horizontal_pixels
        for j = 1:vertical_pixels
            S[i,j,5:8]=e*R2*L*R*S[i,j,5:8]
        end
    end
    return S
    end


function gen_intrinsic_rotation_matrix(pointing=[0,0,0])
    #generate rotation matrix
    #'roll'
    R_x=zeros(4,4)
    R_x[1,1]=1
    R_x[2,2]=1
    R_x[3,3]=cos(pointing[1])
    R_x[4,3]=sin(pointing[1])
    R_x[3,4]=-sin(pointing[1])
    R_x[4,4]=cos(pointing[1])

    #'pitch'
    R_y=zeros(4,4)
    R_y[1,1]=1
    R_y[2,2]=cos(pointing[2])
    R_y[4,2]=-sin(pointing[2])
    R_y[3,3]=1
    R_y[2,4]=sin(pointing[2])
    R_y[4,4]=cos(pointing[2])

    #'yaw'
    R_z=zeros(4,4)
    R_z[1,1]=1
    R_z[2,2]=cos(pointing[3])
    R_z[3,2]=sin(pointing[3])
    R_z[2,3]=-sin(pointing[3])
    R_z[3,3]=cos(pointing[3])
    R_z[4,4]=1

    R=R_z*R_y*R_x
    return R
    end


function gen_final_rotation_matrix(direction=[0,0])
    #turn angle coords into cartesion direction
    d=[sin(direction[1])*cos(direction[2]),
                             sin(direction[1])*sin(direction[2]),
       cos(direction[1])]
    
    #perpendicular vector for rotation
    k=cross([0,0,1],d)
    mag_k = sqrt(transpose(k)*k)
    if mag_k > 0
        k=k/mag_k
        
        #cross product matrix
        K=zeros(3,3)
        K[1,2]=-k[3]
        K[1,3]=k[2]
        K[2,1]=k[3]
        K[2,3]=-k[1]
        K[3,1]=-k[2]
        K[3,2]=k[1]
        
        #make a 3d rotation matrix
        R_3d=Matrix{Float64}(I,3,3)
        R_3d=R_3d + sin(direction[1])*K + (1-cos(direction[1]))*K*K
        
        #make it 4d
        R=zeros(4,4)
        R[1,1]=1
        R[2:4,2:4]=R_3d
        #transpose to give proper handedness
        R=transpose(R)
    else
        println("Final Rotation: Camera velocity already aligned with FIDO z-axis; returning identity matrix")
        R=Matrix{Float64}(I,4,4)
    end
    return R
end


function lorentz_boost_z(beta)
    Lambda=zeros(4,4)
    gamma=1/sqrt(1-beta^2)
    Lambda[1,1]=gamma
    Lambda[4,1]=-beta*gamma
    Lambda[1,4]=-beta*gamma
    Lambda[2,2]=1
    Lambda[3,3]=1
    Lambda[4,4]=gamma
    return Lambda
    end


function calc_lower_metric(position)
    #In this build I'm going with cartesian minkowski space
    g=Matrix{Float64}(I,4,4)
    g[1,1]=-1
    return g
    end

function calc_local_vierbein(position)
    #In this build I'm going with cartesian minkowski space
    return Matrix{Float64}(I,4,4)
    end


function calc_christoffel_udd(position,index)
    #In this build I'm going with cartesian minkowski space
    Christoffel=0
    return Christoffel
    end

function calc_planck(T,nu)
    B_nu=2*h*nu^3/(c^2*(exp(h*nu/(k_B*T))-1))
    return B_nu
end


function get_temp(position)
    #gives the temperature in Kelvin
    return 5778
end

function calc_spectral_emission_coeficient(position_velocity,frequency)
    #j_nu = a_nu*B_nu for thermal emission
    #units are m^-1 with default scale
    j_nu=calc_spectral_absorbtion_coeficient(position_velocity,frequency)*calc_planck(get_temp(position_velocity[1:4]),frequency)/map_scale
    return j_nu
end

function is_fire(position)
    x_in_bottom= -0.02<=position[2]<=0.03
    y_in_left= 0<=position[3]<=0.03
    z_in_block= 0<=position[4]<=0.01

    if x_in_bottom && y_in_left && z_in_block
        pos_in_hole = position[2]>0.01 && position[3]>0.01
        return !(pos_in_hole)
    else
        return false
    end
end

function is_fire2(position)
    if abs(position[2])<=0.01 && abs(position[3])<=0.01 && -0.01<=position[4]<=0.5
        return true
    else
        return false
    end
end

function calc_spectral_absorbtion_coeficient(position_velocity,frequency)
    #units are m^-1 with default scale
    #Set a_nu=760 for a real value from abstract of B.L. Wersborg, L.K. Fox, J.B. Howard 1974
    #not public :(
    if is_fire2(position_velocity[1:4])
        a_nu=50/map_scale
        return a_nu
    else
        return 0
    end
end


function get_source_velocity(position)
    return [1,0,0,0]
end


function calc_ray_derivative(Ray,raylength,colors_freq)
    #calculate derivative at a point
    slope= zeros(raylength)
    for i in 1:4
        #dx/dl=v
        slope[i]=Ray[i+4]
    end
    for i in 5:8
        #geodesic equation
        a=0
        for j in 5:8
            for k in 5:8
                a=a - calc_christoffel_udd(Ray[1:4],CartesianIndex(i-4,j-4,k-4))*Ray[j]*Ray[k]
            end
        end
        slope[i]=a
    end
    for i in 9:raylength
        if isodd(i)
            #calculate the frequency of the ray in the source frame, by nu=E/hbar, E = -p * u
            freq_shift=-transpose(get_source_velocity(Ray[1:4]))*calc_lower_metric(Ray[1:4])*Ray[5:8]
            nu=colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the invariant brightness the ray "will" (hence the - sign) gain between here and the camera
            slope[i]=-calc_spectral_emission_coeficient(Ray[1:8],nu)*exp(-Ray[i+1])/nu^3
        else
            #calculate the frequency of the ray in the source frame, by nu=E/hbar, E = -p * u
            freq_shift=-transpose(get_source_velocity(Ray[1:4]))*calc_lower_metric(Ray[1:4])*Ray[5:8]
            nu=-colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the optical depth which the ray "will" (hence the - sign) pass through in the "future" (+lambda) to get to the camera
            slope[i]=-calc_spectral_absorbtion_coeficient(Ray[1:8],nu)*freq_shift
        end
    end
    return slope
end


function RK4_Step(Ray,raylength,stepsize,colors_freq)
    #Runge-Kutta Method -- Classic
    k1 = calc_ray_derivative(Ray,raylength,colors_freq)
    k2 = calc_ray_derivative(Ray + k1*stepsize/2, raylength,colors_freq)
    k3 = calc_ray_derivative(Ray + k2*stepsize/2, raylength,colors_freq)
    k4 = calc_ray_derivative(Ray + k3*stepsize, raylength,colors_freq)
    return Ray + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
end

function RKDP_Step(Ray,slope_last,raylength::Integer,stepsize::Real,colors_freq,abs_tol,rel_tol,prev_rejected::Bool,
                   max_dt::Real)
    #Runge-Kutta Method -- Dormand-Prince method
    #Calculating k's
    k1 = slope_last
    k2 = calc_ray_derivative(Ray + k1*stepsize*DP_Butcher[2,2], raylength,colors_freq)
    k3 = calc_ray_derivative(Ray + stepsize*(k1*DP_Butcher[3,2] + k2*DP_Butcher[3,3]),raylength,colors_freq)
    k4 = calc_ray_derivative(Ray + stepsize*(k1*DP_Butcher[4,2] + k2*DP_Butcher[4,3] + k3*DP_Butcher[4,4]),
                             raylength,colors_freq)
    k5 = calc_ray_derivative(Ray + stepsize*(k1*DP_Butcher[5,2] + k2*DP_Butcher[5,3] + k3*DP_Butcher[5,4] +
        k4*DP_Butcher[5,5]),raylength, colors_freq)
    k6 = calc_ray_derivative(Ray + stepsize*(k1*DP_Butcher[6,2] + k2*DP_Butcher[6,3] + k3*DP_Butcher[6,4] +
        k4*DP_Butcher[6,5] + k5*DP_Butcher[6,6]),raylength, colors_freq)
    
    next_slope=k1*DP_Butcher[7,2] + k3*DP_Butcher[7,4] + k4*DP_Butcher[7,5] + k5*DP_Butcher[7,6] +
        k6*DP_Butcher[7,7]
    
    k7 = calc_ray_derivative(Ray + stepsize*next_slope,raylength, colors_freq)
    
    #Calculating y's
    y = Ray + stepsize*slope_last
    y_hat = Ray + stepsize*(k1*DP_Butcher[9,2] + k3*DP_Butcher[9,4] + k4*DP_Butcher[9,5] + k5*DP_Butcher[9,6] +
        k6*DP_Butcher[9,7] + k7*DP_Butcher[9,8])

    #Estimate error
    delta_y = y-y_hat
    error=0
    for i in 1:raylength
        tol = abs_tol[i]+max(abs(Ray[i]),abs(y[i]))*rel_tol[i]
        error += (delta_y[i]/tol)^2
    end
    error=sqrt(error/raylength)

    """
    if isnan(stepsize)
        println(true)
    end
    """
    
    #decide next step
    #calc size of next step
    new_stepsize=0.9*stepsize*(1/error)^(1/5)

    if abs(new_stepsize)>max_dt
        new_stepsize=-max_dt
    end
    
    if error == 0 || !isfinite(new_stepsize)
        #error small enough; send updated stats marked as accepted, with same step size (to keep out NaN's)
        return y,next_slope,stepsize,false
    elseif error > 1
        #error too large; send back to manager marked as rejected, with new step size
        if prev_rejected && abs(new_stepsize)>=abs(stepsize)
            #make sure stepsize goes down when the error is too large
            new_stepsize=stepsize/2
        end
        return Ray,slope_last,new_stepsize,true
    else
        #error small enough; send updated stats marked as accepted, with new step size
        #keep stepsize from increasing quickly
        if abs(new_stepsize)>5*abs(stepsize)
            new_stepsize=5*stepsize
        end
        return y,next_slope,new_stepsize,false
    end
end


function Heaviside(x)
    if x >= 0
        return 1
    else
        return 0
    end
end

function Selector_funct(x,y,z)
    return y*(1-Heaviside(x))+z*Heaviside(x)
end

function cie_x(wavelength)
    x=0
    for i in 1:3
        x=x+cie_matrix[1,i]*exp((-1/2)*((wavelength-cie_matrix[2,i])*Selector_funct(
            wavelength-cie_matrix[2,i],cie_matrix[3,i],cie_matrix[4,i]))^2)
    end
    return x
end

function cie_y(wavelength)
    y=0
    for i in 1:2
        y=y+cie_matrix[1,i+3]*exp((-1/2)*((wavelength-cie_matrix[2,i+3])*Selector_funct(
            wavelength-cie_matrix[2,i+3],cie_matrix[3,i+3],cie_matrix[4,i+3]))^2)
    end
    return y
end

function cie_z(wavelength)
    z=0
    for i in 1:2
        z=z+cie_matrix[1,i+5]*exp((-1/2)*((wavelength-cie_matrix[2,i+5])*Selector_funct(
            wavelength-cie_matrix[2,i+5],cie_matrix[3,i+5],cie_matrix[4,i+5]))^2)
    end
    return z
end

function ray_to_I_lambda(ray,colors,colors_freq)
    n_colors=length(colors)
    I_lambdas=zeros(n_colors)
    for i in 1:n_colors
        #I_lambda=I_nu *c/lambda^2, where c=lambda*nu
        #since colors are in nm, this gives the spectral radiance in W sr^-1 m^-2 nm^-1
        I_lambdas[i]=ray[7+2*i]*colors_freq[i]^3*c/colors[i]^2
    end
    return I_lambdas
end

function calc_xyY(ray,colors,colors_freq)
    I_lambdas=ray_to_I_lambda(ray,colors,colors_freq)
    
    I_interpolation=linear_interpolation(colors,I_lambdas)
    
    delXYZ(lambda,p) = I_interpolation(lambda)*[cie_x(lambda),cie_y(lambda),cie_z(lambda)]
    domain = (minimum(colors),maximum(colors))
    prob=IntegralProblem(delXYZ,domain)
    CIEXYZ=solve(prob,HCubatureJL();reltol=1e-3,abstol=1e-3)   
    if CIEXYZ[1] != 0 && CIEXYZ[2] != 0 && CIEXYZ[3] != 0
        CIEXYZ=XYZ{Float64}(CIEXYZ[1],CIEXYZ[2],CIEXYZ[3])
        CIExyY=xyY(CIEXYZ)
    else
        CIExyY=xyY{Float64}(1,1,0)
    end
    
    return CIExyY
end

function alltogethernow(;camera_pos=[0,0,0,-0.05],camera_dir=[0.0,0.0],speed=0.0,
                        camera_point=[0.0,0.0,0.0],x_pix=40,stepsize=0.0001,n_steps=1000,colors=[400,550,700],
                        returnrays=false)
    #initialize rays
    ray_matrix=initialize_camera(camera_pos,direction=camera_dir,beta=speed,pointing=camera_point,
                                 horizontal_pixels=x_pix,colors=colors)
    #initialize xyY pixels
    xyY_img=Array{xyY{Float64}}(undef,size(ray_matrix,1),size(ray_matrix,2))
    
    raylength=8+2*length(colors)
    f(x)=c*1e9/x

    #initialize colors_freq
    colors_freq=f.(colors)
    
    for i in 1:size(ray_matrix,1)
        for j in 1:size(ray_matrix,2)
            ray=ray_matrix[i,j,:]
            #integrate ray
            for k in 1:n_steps
                ray=RK4_Step(ray,raylength,-stepsize,colors_freq)               
            end

            if returnrays
                ray_matrix[i,j,:]=ray
            end
            #calculate pixel value
            xyY_img[i,j]=calc_xyY(ray,colors,colors_freq)
        end
    end

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

function alltogethernow_RKDP(;camera_pos=[0,0,0,0],camera_dir=[0.0,0.0],speed=0.0,
                        camera_point=[0.0,0.0,0.0],x_pix=40,dt=1e-3,max_steps=1e3,colors=[400,550,700],
                        returnrays=false,tolerance=1e-8,max_dt=1e-1,fov_hor=85::Real, fov_vert=60::Real)
    #initialize rays
    ray_matrix=initialize_camera(camera_pos,direction=camera_dir,beta=speed,pointing=camera_point,
                                 horizontal_pixels=x_pix,colors=colors,fov_hor=fov_hor, fov_vert=fov_vert)
    #initialize xyY pixels
    xyY_img=Array{xyY{Float64}}(undef,size(ray_matrix,1),size(ray_matrix,2))

    n_colors=length(colors)
    raylength=8+2*n_colors
    f(x)=c*1e9/x

    #initialize colors_freq
    colors_freq=f.(colors)

    #make tolerance arrays
    abs_tol=fill(1e-4,length(ray_matrix[1,1,:]))
    rel_tol=fill(1e-4,length(ray_matrix[1,1,:]))

    #flip to back in time
    dt=-dt

    for i in 1:size(ray_matrix,1)
        for j in 1:size(ray_matrix,2)
            inner_dt=dt
            ray=ray_matrix[i,j,:]
            #integrate ray
            raytrace=true
            shared_slope=calc_ray_derivative(ray,raylength,colors_freq)
            rejected=false
            step_count=0
            
            while raytrace
                dt_old=inner_dt
                ray,shared_slope,inner_dt,rejected=RKDP_Step(ray,shared_slope,raylength,inner_dt,colors_freq,
                                                       abs_tol, rel_tol,rejected,max_dt)
                step_count+=1
                #my current termination condition
                min_optical_depth=minimum(ray[9:2:end])
                dist_from_origin=sqrt(ray[2]^2 + ray[3]^2 + ray[4]^2)
                if step_count>=max_steps || min_optical_depth>=-log(tolerance) || dist_from_origin > 0.7
                    #println(step_count," ",dt_old," ",inner_dt)
                    raytrace=false
                end
            end
            
            if returnrays
                ray_matrix[i,j,:]=ray
            end
            #calculate pixel value
            xyY_img[i,j]=calc_xyY(ray,colors,colors_freq)
        end
    end

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
