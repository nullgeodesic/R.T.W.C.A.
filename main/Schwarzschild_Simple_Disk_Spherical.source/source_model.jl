"""
Author: Levi Malmström
"""

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
@inline function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    terminate = false
    if step_count >= max_steps || ray[2] > 100 || ray[2] <= (10*rel_tol[2] + 2)*r_s || ray[4] > 8π
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
