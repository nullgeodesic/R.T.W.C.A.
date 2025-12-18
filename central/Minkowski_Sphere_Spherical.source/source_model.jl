"""
v0.3.5
December 12 2025
Author: Levi Malmström
"""

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
    @views j_nu=calc_spectral_absorbtion_coeficient(position_velocity,frequency)*calc_planck(get_temp(position_velocity[1:4]),frequency)
    return j_nu
end

function calc_spectral_absorbtion_coeficient(position_velocity,frequency)
    #units are m^-1 with default scale
    @views if is_fire(position_velocity[1:4])
        a_nu=1/map_scale
        return a_nu
    else
        return 0
    end
end

function is_fire(position)
    #sphere of radius 1 centered on origin
    if position[2] <=1
        return true
    else
        return false
    end
end

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


function get_source_velocity(position)
    return [1,0,0,0]
end

function calc_terminate(ray,dt,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale, max_steps,step_count)
    @views if step_count >= max_steps || minimum(ray[9:2:end]) >=
        -log(abs_tol[10]) || ray[2] > 100
        return true
    else
        return false
    end
end
