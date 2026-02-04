"""
Author: Levi Malmström
"""

#CONSTANTS
#emission scale factor (to scale the rays at the end)
#2h/c^2, but scaled to ν measured in 1/nm instead of Hz
const emission_scale = 2*h*c*1e27
#h/k_B in nm/K, instead of the normal s/K, because ν is in units of nm^-1, not Hz
const ν_factor = h*c*1e9/k_B



"""
Planck function (ν in nm^-1).
"""
@inline function calc_planck(T,ν)
    B_ν = emission_scale*ν^3/(exp(ν_factor*ν/T)-1)
    return B_ν
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
    if position[2] <=1
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
    max_step_size = min(inv(abs(inv(max_dt_scale/(abs(ray[6])+no_div_zero)) +
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
    if step_count >= max_steps || ray[2] > 100 || ray[4] > 8π
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


function skybox_handling!(ray,raylength,colors,colors_freq,beamsize)
    #no skybox
    return nothing
end
