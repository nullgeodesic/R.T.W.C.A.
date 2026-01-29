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
