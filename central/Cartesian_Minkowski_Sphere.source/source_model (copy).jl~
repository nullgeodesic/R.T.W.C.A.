"""
v0.3.2
October 3 2025
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
    @views j_nu=calc_spectral_absorbtion_coeficient(position_velocity,frequency)*calc_planck(get_temp(position_velocity[1:4]),frequency)/map_scale
    return j_nu
end

function calc_spectral_absorbtion_coeficient(position_velocity,frequency)
    #units are m^-1 with default scale
    #Set a_nu=760 for a real value from abstract of B.L. Wersborg, L.K. Fox, J.B. Howard 1974
    #not public :(
    @views if is_fire3(position_velocity[1:4])
        a_nu=50/map_scale
        return a_nu
    else
        return 0
    end
end

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

function pad_max_dt2(position,max_dt_scale)
    if abs(position[2])<=0.02 && abs(position[3])<=0.02 && -0.01<=position[4]<=0.11
        return max_dt_scale
    else
        dist=sqrt((abs(position[2])-0.02)^2 + (abs(position[3])-0.02)^2 + min(abs(position[4] + 0.01),
                                                                              abs(position[4]-0.11))^2)
        return max_dt_scale*(1 + dist)
    end
end

function pad_max_dt3(position,max_dt_scale)
    return max_dt_scale*(1+abs(position[2])+abs(position[3])+abs(position[4]))
end

function get_source_velocity(position)
    return [1,0,0,0]
end
