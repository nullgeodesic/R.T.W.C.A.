"""
Author: Levi Malmström
"""
#CONSTANTS
#Skybox 1
const skybox1 = load("Textures/Real_Skies/NASA_Panorama_1024.png")
#height of skybox1 in pixels
const skybox1_pix_height = size(skybox1,2) ÷ 4
#Skybox 2
const skybox2 = load("Textures/Nebulas/Red_Nebula.png")
#height of skybox2 in pixels
const skybox2_pix_height = size(skybox2,2) ÷ 4

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
@inline function get_temp(ray)
    if ray[2] >0
        return 5778
    else
        return 4000
    end
end


"""
Calculates the spectral emission coeficient for a BB radiator.
"""
@inline function calc_spectral_emission_coeficient(ray,frequency)
    #j_nu = a_nu*B_nu for thermal emission
    #units are m^-1 with default scale
    j_nu = calc_spectral_absorbtion_coeficient(ray,frequency)*calc_planck(get_temp(ray),frequency)
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
    """
    #wedge shaped disk from 6 <= r <= 12, 15π/32 <= θ <= 17π/32
    if 6 <= r_of_l(position[2]) <= 12 && 15π/32 <= position[3] <= 17π/32
        return true
    else
        return false
    end
    """
    return false
end


"""
Keeps the integrator from going too fast.
"""
@inline function pad_max_dt(ray,max_dt_scale)
    #dλ = 1/(|dλ1|^-1 + |dλ2|^-1 + |dλ3|^-1)
    #dλ1 = ϵ/(|v| + δ)
    #dλ2 = ϵ min(θ,π-θ)/(|ω| + δ)
    #dλ3 = ϵ/(|ψ| + δ)
    max_step_size = min(inv(abs(inv(max_dt_scale/(abs(ray[6]) + no_div_zero)) +
        no_div_zero) + abs(inv(max_dt_scale*min(ray[3],π-ray[3])/(abs(ray[7]) +
        no_div_zero) + no_div_zero)) +
        abs(inv(max_dt_scale/(abs(ray[8])+no_div_zero) + no_div_zero)) + no_div_zero),0.5)
    max_step_size = max(max_step_size,1e-12)
    return max_step_size
end


"""
Velocity of the material at a position (mutating).
"""
@inline function get_source_velocity!(position,source_vel_buffer)
    #FIDO
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
    return [1.0,0.0,0.0,0.0]
end


"""
Whether to stop integrating the ray.
"""
@inline function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    terminate = false
    if step_count >= max_steps || r_of_l(ray[2]) > 100 || ray[4] > 8π
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


function skybox_I_ν_calc(i::Int,j::Int,f::Real,beamsize::Real,skybox_num::Int)
    #Note: pretending sRGB/Rec. 709 is the same as Rec. 2020 for this crude RGB to spectrum scheme
    #Also, j is flipped in the skybox1 to j_new =  3*skybox1_pix_height + 1 - j
    if skybox_num == 1
        R = Float64(skybox1[3*skybox1_pix_height + 1 - j,i].r)
        #λ = 630 nm; f = 0.001587 nm^-1
        G = 0.5*Float64(skybox1[3*skybox1_pix_height + 1 - j,i].g)
        #λ = 532 nm; f = 0.001879 nm^-1
        B = 0.32*Float64(skybox1[3*skybox1_pix_height + 1 - j,i].b)
        #λ = 467 nm; f = 0.002141 nm^-1
    elseif skybox_num == 2
        R = Float64(skybox2[3*skybox2_pix_height + 1 - j,i].r)
        #λ = 630 nm; f = 0.001587 nm^-1
        G = 0.5*Float64(skybox2[3*skybox2_pix_height + 1 - j,i].g)
        #λ = 532 nm; f = 0.001879 nm^-1
        B = 0.32*Float64(skybox2[3*skybox2_pix_height + 1 - j,i].b)
        #λ = 467 nm; f = 0.002141 nm^-1
    end
    #line width ~ 10 nm
    return 1e-12 * (R*Normal_Gauss(f,2e-5,0.001587) + G*Normal_Gauss(f,2e-5,0.001879) +
        B*Normal_Gauss(f,2e-5,0.002141))
end


function skybox_handling!(ray,raylength,colors,colors_freq,beamsize)
    if r_of_l(ray[2]) > 50 && ray[2] > 0
        skybox_num = 1
        position4 = ray[1:4]
        velocity4 = ray[5:8]
        g = calc_lower_metric(position4)
        freq_shift = -transpose([1,0,0,0])*g*velocity4
        #transform to flat cartesian coordinates
        e = calc_vierbein(position4)
        R = R_sph_to_cart(position4)
        v = R*e*velocity4
        
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
        i = min(floor(Int64,skybox1_pix_height*(u + 1)/2) + 1,skybox1_pix_height)
        j = min(floor(Int64,skybox1_pix_height*(v + 1)/2) + 1,skybox1_pix_height)
        #shift to correct face
        if face == 1
            j += skybox1_pix_height
        elseif face == 2
            i += skybox1_pix_height
            j += skybox1_pix_height
        elseif face == 3
            i += 2*skybox1_pix_height
            j += skybox1_pix_height
        elseif face == 4
            i += 3*skybox1_pix_height
            j += skybox1_pix_height
        elseif face == 5
            i += skybox1_pix_height
        else
            i += skybox1_pix_height
            j += 2*skybox1_pix_height
        end
            
        for k in 9:2:raylength
            f = colors_freq[ceil(Int,(k-8)/2)]*freq_shift
            ray[k] += skybox_I_ν_calc(i,j,f,beamsize,skybox_num)*exp(-ray[k+1])/f^3
        end


    elseif r_of_l(ray[2]) > 50 && ray[2] < 0
        skybox_num = 2
        position4 = ray[1:4]
        velocity4 = ray[5:8]
        g = calc_lower_metric(position4)
        freq_shift = -transpose([1,0,0,0])*g*velocity4
        #transform to flat cartesian coordinates
        e = calc_vierbein(position4)
        R = R_sph_to_cart(position4)
        v = R*e*velocity4
        
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
        i = min(floor(Int64,skybox2_pix_height*(u + 1)/2) + 1,skybox2_pix_height)
        j = min(floor(Int64,skybox2_pix_height*(v + 1)/2) + 1,skybox2_pix_height)
        #shift to correct face
        if face == 1
            j += skybox2_pix_height
        elseif face == 2
            i += skybox2_pix_height
            j += skybox2_pix_height
        elseif face == 3
            i += 2*skybox2_pix_height
            j += skybox2_pix_height
        elseif face == 4
            i += 3*skybox2_pix_height
            j += skybox2_pix_height
        elseif face == 5
            i += skybox2_pix_height
        else
            i += skybox2_pix_height
            j += 2*skybox2_pix_height
        end
            
        for k in 9:2:raylength
            f = colors_freq[ceil(Int,(k-8)/2)]*freq_shift
            ray[k] += skybox_I_ν_calc(i,j,f,beamsize,skybox_num)*exp(-ray[k+1])/f^3
        end
    end
    return nothing
end
