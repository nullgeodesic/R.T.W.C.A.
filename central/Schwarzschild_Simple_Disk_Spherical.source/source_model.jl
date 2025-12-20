"""
v0.3.7
December 19 2025
Author: Levi Malmström
"""

"""
Planck function.
"""
function calc_planck(T,nu)
    B_nu=2*h*nu^3/(c^2*(exp(h*nu/(k_B*T))-1))
    return B_nu
end


"""
Gives the temperature at a position in Kelvin.
"""
function get_temp(position)
    return 5778
end


"""
Calculates the spectral emission coeficient for a BB radiator.
"""
function calc_spectral_emission_coeficient(pos_vel,frequency)
    #j_nu = a_nu*B_nu for thermal emission
    #units are m^-1 with default scale
    @views j_nu=calc_spectral_absorbtion_coeficient(pos_vel,frequency)*calc_planck(get_temp(pos_vel),frequency)
    return j_nu
end


"""
Calculates the spectral absorbption coeficient.
"""
function calc_spectral_absorbtion_coeficient(pos_vel,frequency)
    #units are 1/M
    @views if is_fire(pos_vel)
        a_nu=1/map_scale
        return a_nu
    else
        return 0
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
function pad_max_dt(pos_vel,max_dt_scale)
    #dλ = 1/(|dλ1|^-1 + |dλ2|^-1 + |dλ3|^-1)
    #dλ1 = ϵ/(|v| + δ)
    #dλ2 = ϵ min(θ,π-θ)/(|ω| + δ)
    #dλ3 = ϵ/(|ψ| + δ)
    max_step_size = min(inv(abs(inv(max_dt_scale/(abs(pos_vel[6])+no_div_zero)) +
        no_div_zero) + abs(inv(max_dt_scale*min(pos_vel[3],π-pos_vel[3])/(abs(pos_vel[7]) +
        no_div_zero) + no_div_zero)) +
        abs(inv(max_dt_scale/(abs(pos_vel[8])+no_div_zero) + no_div_zero)) + no_div_zero),0.5)
    max_step_size = max(max_step_size,1e-12)
    return max_step_size
end


"""
Velocity of the material at a position (mutating).
"""
function get_source_velocity!(position,source_vel_buffer::Vector)
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
