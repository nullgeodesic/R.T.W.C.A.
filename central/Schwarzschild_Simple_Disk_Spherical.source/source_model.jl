"""
v0.4.0
December 24 2025
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
function calc_planck(T,ν)
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
function get_temp(ray)
    return 5778
end


"""
Calculates the spectral emission coeficient for a BB radiator.
"""
function calc_spectral_emission_coeficient(ray,frequency)
    #j_nu = a_nu*B_nu for thermal emission
    #units are m^-1 with default scale
    j_nu=calc_spectral_absorbtion_coeficient(ray,frequency)*calc_planck(get_temp(ray),frequency)
    return j_nu
end


"""
Calculates the spectral absorbption coeficient.
"""
function calc_spectral_absorbtion_coeficient(ray,frequency)
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
function is_fire(position)
    #wedge shaped disk from 6 <= r <= 12, 15π/32 <= θ <= 17π/32
    if 6 <= position[2] <= 12 && 15π/32 <= position[3] <= 17π/32
        return true
    else
        return false
    end
end


"""
Keeps the integrator from going too fast.
"""
function pad_max_dt(ray,max_dt_scale)
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
function get_source_velocity!(position,source_vel_buffer)
    #M=1
    #Ω = sqrt(M/r^3)
    #ut = (1 - 3M/r)^(-1/2)
    #uα = ut*[1,0,0,Ω]
    if is_fire(position)
        #circular orbit
        source_vel_buffer[1] = inv(sqrt(1 - 3*M/position[2]))
        source_vel_buffer[2] = 0.0
        source_vel_buffer[3] = 0.0
        source_vel_buffer[4] = inv(sqrt(1 - 3*M/position[2]))*sqrt(M/position[2]^3)
        return nothing
    else
        #FIDO
        source_vel_buffer[1] = inv(sqrt(1 - 2*M/position[2]))
        source_vel_buffer[2] = 0.0
        source_vel_buffer[3] = 0.0
        source_vel_buffer[4] = 0.0
        return nothing
    end

end


"""
Velocity of the material at a position (non-mutating).
"""
function get_source_velocity(position)
    #M=1
    #Ω = sqrt(M/r^3)
    #ut = (1 - 3M/r)^(-1/2)
    #uα = ut*[1,0,0,Ω]
    if is_fire(position)
        return [inv(sqrt(1 - 3*M/position[2])),0,0,inv(sqrt(1 - 3*M/position[2]))*sqrt(M/position[2]^3)]
    else
        return [inv(sqrt(1 - 2*M/position[2])),0,0,0]
    end

end


"""
Whether to stop integrating the ray.
"""
function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    @views if step_count >= max_steps || minimum(ray[9:2:end]) >=
        -log(abs_tol[10]) || ray[2] > 100 || ray[2] <= (10*rel_tol[2] + 2)*r_s || ray[4] > 8π
        return true
    else
        return false
    end
end
