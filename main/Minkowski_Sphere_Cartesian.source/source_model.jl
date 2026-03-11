"""
Author: Levi Malmström
"""

#CONSTANTS


function load_textures()
    #Skybox 1
    skybox1 = load("Textures/Test_Cubemaps/Cubemap_Sky_22-512x512.png")
    #height of skybox1 in pixels
    skybox1_pix_height = size(skybox1,2) ÷ 4
    #Skybox 2
    skybox2 = 1
    #height of skybox2 in pixels
    skybox2_pix_height = 1
    return skybox1,skybox1_pix_height,skybox2,skybox2_pix_height
end


"""
Planck function (f in nm^-1).
"""
@inline function calc_planck(T,f)
    B_f = emission_scale*f^3/(exp(f_factor*f/T)-1)
    return B_f
end

"""
Planck function (ν in Hz).
"""
function calc_planck_Hz(T,ν)
    B_ν = 2*h*ν^3/(c^2*(exp(h*ν/(k_B*T))-1))
    return B_ν
end


"""
Gives the temperature at a position in Kelvin.
"""
@inline function get_temp(position)
    #gives the temperature in Kelvin
    return 5778
end


"""
Calculates the spectral emission coeficient for a BB radiator.
"""
@inline function calc_spectral_emission_coeficient(ray,frequency)
    #j_nu = a_nu*B_nu for thermal emission
    #units are m^-1 with default scale
    j_nu=calc_spectral_absorbtion_coeficient(ray,frequency)*calc_planck(get_temp(ray),frequency)
    return j_nu
end


"""
Calculates the spectral absorbption coeficient.
"""
@inline function calc_spectral_absorbtion_coeficient(ray,frequency)
    #units are 1/M
    a_nu = 1/map_scale
    if is_fire(ray)
        return a_nu
    else
        return 0.0
    end
end



"""
Determines if there is emmiting material at a location.
"""
@inline function is_fire(position)
    #sphere of radius 1 centered on origin
    if sqrt(position[2]^2 + position[3]^2 + position[4]^2) <= 1
        return true
    else
        return false
    end
end


"""
function is_fire(position)
    #lopsided L
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
    #rectangular prism
    if abs(position[2])<=0.01 && abs(position[3])<=0.01 && 0<=position[4]<=0.1
        return true
    else
        return false
    end
end

function is_fire3(position)
    #sphere
    if sqrt(position[2]^2 + position[3]^2 + position[4]^2) <= 0.02
        return true
    else
        return false
    end
end
"""


"""
Keeps the integrator from going too fast.
"""
@inline function pad_max_dt(position,max_dt_scale)
    return max_dt_scale*(1+abs(position[2])+abs(position[3])+abs(position[4]))
end


"""
Velocity of the material at a position (mutating).
"""
@inline function get_source_velocity!(position,source_vel_buffer)
    source_vel_buffer[1] = 1.0
    source_vel_buffer[2] = 0.0
    source_vel_buffer[3] = 0.0
    source_vel_buffer[4] = 0.0
    return nothing
end


"""
Velocity of the material at a position (non-mutating).
"""
function get_source_velocity(position)
    return [1,0,0,0]
end


"""
Whether to stop integrating the ray.
"""
@inline function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    terminate = false
    if step_count >= max_steps || sqrt(ray[2]^2 + ray[3]^2 + ray[4]^2) > 100
        terminate = true
    else
        terminate = true
        for i in 10:2:raylength
            if ray[i] < -log(abs_tol[10])
                terminate = false
            end
        end
    end

    return terminate
end


function skybox_I_ν_calc(i::Integer,j::Integer,f::Real,skybox_num::Integer,skybox1,skybox1_pix_height,skybox2,
                         skybox2_pix_height)
    #Note: pretending sRGB/Rec. 709 is the same as Rec. 2020 for this crude RGB to spectrum scheme
    #Also, j is flipped in the skybox to j_new =  3*skybox_pix_height + 1 - j
    R = Float64(skybox1[3*skybox1_pix_height + 1 - j,i].r)
    #λ = 630 nm; f = 0.001587 nm^-1
    G = 0.5*Float64(skybox1[3*skybox1_pix_height + 1 - j,i].g)
    #λ = 532 nm; f = 0.001879 nm^-1
    B = 0.32*Float64(skybox1[3*skybox1_pix_height + 1 - j,i].b)
    #λ = 467 nm; f = 0.002141 nm^-1
    #line width ~ 10 nm
    return 1e-12 * (R*Normal_Gauss(f,2e-5,0.001587) + G*Normal_Gauss(f,2e-5,0.001879) +
        B*Normal_Gauss(f,2e-5,0.002141))
end


function calc_skybox_coord(position4,velocity4,skybox_pix_height::Int)
    #space already flat and cartesian
    """
    #transform to flat cartesian coordinates
    e = calc_vierbein(position4)
    R = R_sph_to_cart(position4)
    v = R*e*velocity4
    """
    v = velocity4
    
    #start by getting a normalized version of the 3 velocity (could get by dividing by γ, but this way
    #should be robust against floating point errors)
    v3 = v[2:4]/sqrt(v[2]^2 + v[3]^2 + v[4]^2)

    x = v3[1]
    y = v3[2]
    z = v3[3]
    #figure out which face of the skybox the ray is going to hit
    if abs(x) >= abs(y) && abs(x) >= abs(z)
        if x > 0
            face = 1
        else
            face = 3
        end
    elseif abs(y) >= abs(x) && abs(y) >= abs(z)
        if y > 0
            face = 4
        else
            face = 2
        end
    else
        if z > 0
            face = 6
        else
            face = 5
        end
    end

    #project ray onto skybox
    if face == 1
        #w = -y
        #h = z
        #maj = abs(x)
        #u = w/maj
        #v = h/maj
        u = -y/abs(x)
        v = z/abs(x)
    elseif face == 2
        #w = -x
        #h = z
        #maj = abs(y)
        #u = w/maj
        #v = h/maj
        u = -x/abs(y)
        v = z/abs(y)
    elseif face == 3
        #w = y
        #h = z
        #maj = abs(x)
        u = y/abs(x)
        v = z/abs(x)
    elseif face == 4
        #w = x
        #h = z
        #maj = abs(y)
        u = x/abs(y)
        v = z/abs(y)
    elseif face == 5
        #w = -x
        #h = -y
        #maj = abs(z)
        u = -x/abs(z)
        v = -y/abs(z)
    elseif face == 6
        #w = -x
        #h = y
        #maj = abs(z)
        u = -x/abs(z)
        v = y/abs(z)
    end

    #convert from uv into ij (pixels)
    i = min(floor(Int64,skybox_pix_height*(u + 1)/2) + 1,skybox_pix_height)
    j = min(floor(Int64,skybox_pix_height*(v + 1)/2) + 1,skybox_pix_height)
    #note: i,j are local to the face
    """
    if face == 1
        j += skybox_pix_height
    elseif face == 2
        i += skybox_pix_height
        j += skybox_pix_height
    elseif face == 3
        i += 2*skybox_pix_height
        j += skybox_pix_height
    elseif face == 4
        i += 3*skybox_pix_height
        j += skybox_pix_height
    elseif face == 5
        i += skybox_pix_height
    else
        i += skybox_pix_height
        j += 2*skybox_pix_height
    end
    """
    return i,j,u,v,face
end


function calc_skybox_face(i::Integer,j::Integer,skybox_pix_height::Integer)
    if 1 <= i <= 2*skybox_pix_height
        if 1 <= i <= skybox_pix_height
            if skybox_pix_height < j <= 2*skybox_pix_height
                face = 1
            else
                face = 0
            end
        elseif skybox_pix_height < i <= 2*skybox_pix_height
            if 1 <= j <= skybox_pix_height
                face = 5
            elseif skybox_pix_height < j <= 2*skybox_pix_height
                face = 2
            else
                face = 6
            end
        else
            face = 0
        end
    elseif 2*skybox_pix_height < i <= 3*skybox_pix_height && skybox_pix_height < j <= 2*skybox_pix_height
        face = 3
    elseif 3*skybox_pix_height < i <= 4*skybox_pix_height && skybox_pix_height < j <= 2*skybox_pix_height
        face = 4
    else
        face = 0
    end
    return face
end



function calc_skybox_I_ν_beam!(ray,freq_shift::Real,i::Integer,j::Integer,u::Real,v::Real,face::Integer,
                               raylength::Integer,skybox_pix_height::Integer, skybox_num::Integer, colors_freq,
                               skybox1,skybox1_pix_height,skybox2,skybox2_pix_height)
    #calculate the ray beam shape
    
    δ_prime = ray[raylength]*freq_shift
    
    δ_major = δ_prime*(sqrt(ray[raylength-5]^2 + ray[raylength-4]^2)
                       + sqrt(ray[raylength-3]^2 + ray[raylength-2]^2))
    δ_minor = δ_prime*(sqrt(ray[raylength-5]^2 + ray[raylength-4]^2)
                       - sqrt(ray[raylength-3]^2 + ray[raylength-2]^2))
    μ = ray[raylength-6] + angle(complex(ray[raylength-5],ray[raylength-4])*
        complex(ray[raylength-3],ray[raylength-2]))/2

    #approximate angular size of center-face skybox pixels
    sky_beamscale = atan(2/(skybox_pix_height + 1))/2
    #Ideally, I would scale the skybox angular size by u,v position, but the method in quotes was creating artifacts
    #So I simplified it to one beamscale for now.
    """
    #skybox distortions
    θ = atan(sqrt(u^2 + v^2))
    if u > 0
        ϕ = atan(v/u)
    elseif u < 0 && v >= 0
        ϕ = atan(v/u) + π
    elseif u < 0 && v < 0
        ϕ = atan(v/u) - π
    elseif u == 0 && v > 0
        ϕ = π/2
    elseif u == 0 && v < 0
        ϕ = -π/2
    else
        ϕ = 0
    end

    Δθ_r = sky_beamscale*cos(θ)^2
    Δθ_t = sky_beamscale*abs(cos(θ))
    #Δθ_x and Δθ_y the angular scales of pixels in the x and y direction near the central pixel
    Δθ_x = Δθ_r*abs(cos(ϕ)) + Δθ_t*abs(sin(ϕ))
    Δθ_y = Δθ_r*abs(sin(ϕ)) + Δθ_t*abs(cos(ϕ))
    #figure out beam shape in terms of pixels
    σ_x_new,σ_y_new,μ_new = distort_gauss(δ_major/(δ_minor + no_div_zero),1,μ,δ_prime/(Δθ_x + no_div_zero),
                                          δ_prime/(Δθ_y + no_div_zero))
    """

    σ_x_new,σ_y_new,μ_new = distort_gauss(δ_major,δ_minor,μ,inv(sky_beamscale + no_div_zero),
                                          inv(sky_beamscale + no_div_zero))

    #figure out how far out to collect pixels
    h = min(ceil(Integer,max(σ_x_new,σ_y_new)),10)
    
    #collect weighted spatial average I_ν for each frequency over 4*h^2 pixels
    I_ν = zeros(length(colors_freq))
    for k in -h:h
        for l in -h:h
            #i,j,face for the pixel
            i_wrapped,j_wrapped,face_wrapped = cube_wrap(i,j,face,k,l,skybox_pix_height)
            #absolute cube map coordinates
            i_abs,j_abs = cubemap_coord_calc(i_wrapped,j_wrapped,face_wrapped,skybox_pix_height)
            for m in 9:2:(raylength - 11)
                f = colors_freq[ceil(Int,(m-8)/2)]*freq_shift
                I_ν[ceil(Int,(m-8)/2)] += D2_Gauss(k,l,σ_x = σ_x_new,σ_y = σ_y_new,θ = μ_new)*
                    skybox_I_ν_calc(i_abs,j_abs,f,skybox_num,skybox1,skybox1_pix_height,skybox2,
                                    skybox2_pix_height)*exp(-ray[m+1])/f^3
            end
        end
    end

    for a in 9:2:(raylength - 11)
        ray[a] = I_ν[ceil(Int,(a-8)/2)]
    end

    return nothing
end


function skybox_handling!(ray,raylength::Integer,colors,colors_freq,n_bundle_param,skybox1,skybox1_pix_height,
                          skybox2,skybox2_pix_height)
    skybox_num = 1
    position4 = ray[1:4]
    velocity4 = ray[5:8]
    g = calc_lower_metric(position4)
    freq_shift = -transpose([1,0,0,0])*g*velocity4
    #figure out where on the skybox the ray is
    i,j,u,v,face = calc_skybox_coord(position4,velocity4,skybox1_pix_height)
    #add skybox pixels to the ray
    #use simple vs. bundle rays
    if n_bundle_param == 0
        #absolute cube map coordinates
        i_abs,j_abs = cubemap_coord_calc(i,j,face,skybox1_pix_height)
        for k in 9:2:raylength
            f = colors_freq[ceil(Int,(k-8)/2)]*freq_shift
            ray[k] += skybox_I_ν_calc(i_abs,j_abs,f,skybox_num,skybox1,skybox1_pix_height,skybox2,
                                      skybox2_pix_height)*exp(-ray[k+1])/f^3
        end
    else
        calc_skybox_I_ν_beam!(ray,freq_shift,i,j,u,v,face,
                              raylength,skybox1_pix_height,skybox_num,colors_freq,skybox1,skybox1_pix_height,
                              skybox2,skybox2_pix_height)
    end
    return nothing
end
