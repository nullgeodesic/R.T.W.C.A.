"""
Author: Levi Malmström
"""


function load_textures()
    #No skybox
    return 1,1,1,1
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
@inline function get_temp(ray)
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
    #units are 1/map_scale
    a_nu = 0.1f0/Float32(map_scale)
    if is_fire(ray)
        return a_nu
    else
        return 0.0f0
    end
end


"""
Determines if there is emmiting material at a location.
"""
@inline function is_fire(position)
    #wedge shaped disk in r_ISCO_prograde <= r <= 12, 31π/64 <= θ <= 33π/64
    if r_ISCO_prograde <= position[2] <= 12 && 31π/64 <= position[3] <= 33π/64
        return true
    else
        return false
    end
end


"""
Keeps the integrator from going too fast.
"""
@inline function pad_max_dt(ray,max_dt_scale)
    #dλ = 1/(|dλ1|^-1 + |dλ2|^-1 + |dλ3|^-1)
    #dλ1 = ϵ/(|v| + δ)
    #dλ2 = ϵ min(θ,π-θ)/(|ω| + δ)
    #dλ3 = ϵ/(|ψ| + δ)
    #slow down near poles
    if ray[3] < π/8 || ray[3] > 7π/8
        ϵ = min(Float32(max_dt_scale),1f-1)
    else
        ϵ = Float32(max_dt_scale)
    end
    max_step_size = min(inv(abs(inv(ϵ/(abs(ray[6])+no_div_zero)) +
        no_div_zero) + abs(inv(ϵ*min(ray[3],π-ray[3])/(abs(ray[7]) +
        no_div_zero) + no_div_zero)) +
        abs(inv(ϵ/(abs(ray[8])+no_div_zero) + no_div_zero)) + no_div_zero),0.5f0)
    max_step_size = max(max_step_size,1f-12)
    return max_step_size
end


"""
Velocity of the material at a position (mutating).
"""
@inline function get_source_velocity!(position,source_vel_buffer)
    if is_fire(position)
        #prograde equatorial circular orbit
        Ω = sqrt(M)/(position[2]^(3/2) + α*M)
        D = position[2]^(3/2) - 3*M*sqrt(position[2]) + 2*α*sqrt(M)
        ut = (position[2]^(3/2) + α*sqrt(M))/(position[2]^(3/4)*sqrt(D))
        #uα = ut*[1,0,0,Ω]
        source_vel_buffer[1] = ut
        source_vel_buffer[2] = 0.0f0
        source_vel_buffer[3] = 0.0f0
        source_vel_buffer[4] = ut*Ω
        return nothing
    else
        #FIDO
        A = sqrt(position[2]^2 + (α*cos(position[3]))^2)
        B = position[2]^2 + α^2 - 2*M*position[2]
        C = (position[2]^2 + α^2)^2 - B*(α*sin(position[3]))^2
        #ut = 1/lapse rate
        ut = sqrt(C)/(A*sqrt(B))
        Ω = 2*α*M*position[2]/C
        #uα = ut*[1,0,0,Ω]
        source_vel_buffer[1] = ut
        source_vel_buffer[2] = 0.0f0
        source_vel_buffer[3] = 0.0f0
        source_vel_buffer[4] = ut*Ω
        return nothing
    end

end


"""
Velocity of the material at a position (non-mutating).
"""
function get_source_velocity(position)
    if is_fire(position)
        #prograde equatorial circular orbit
        Ω = sqrt(M)/(position[2]^(3/2) + α*M)
        D = position[2]^(3/2) - 3*M*sqrt(position[2]) + 2*α*sqrt(M)
        ut = (position[2]^(3/2) + α*sqrt(M))/(position[2]^(3/4)*sqrt(D))
        #uα = ut*[1,0,0,Ω]
        return [ut,0,0,ut*Ω]
    else
        #FIDO
        A = sqrt(position[2]^2 + (α*cos(position[3]))^2)
        B = position[2]^2 + α^2 - 2*M*position[2]
        C = (position[2]^2 + α^2)^2 - B*(α*sin(position[3]))^2
        #ut = 1/lapse rate
        ut = sqrt(C)/(A*sqrt(B))
        Ω = 2*α*M*position[2]/C
        #uα = ut*[1,0,0,Ω]
        return [ut,0,0,ut*Ω]
    end

end


"""
Whether to stop integrating the ray.
"""
@inline function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    terminate = false
    if step_count >= max_steps || ray[2] > 100 || ray[2] <= (1 + 10*rel_tol[2])*r_ev || ray[4] > 8π
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


function skybox_handling!(ray,raylength::Integer,colors,colors_freq,n_bundle_param,skybox1,skybox1_pix_height,
                          skybox2,skybox2_pix_height)
    #No skybox
    return nothing
end
