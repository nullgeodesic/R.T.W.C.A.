"""
Author: Levi Malmström
"""

#CONSTANTS
#How many meters corespond to one unit of the map
#map_scale = Mass of M87* in meters
const map_scale = 9.62f12
#Prefactor for synchrotron emissivity
const synch_em = Float32(7.694e-37 * map_scale * 2.8e7)
#Prefactor fo synchrotron absorptivity
const synch_ab = Float32(2.8173e-24 * map_scale)
#Jet electron density scale (10^7 electrons/m^3 = 10 electrons/cm^3)
const ne_0 = 1f7

#Power Law Synchrotron:
synch_integrand(μ,p) = (1 + (p[2] - 1)*μ^2)^(-p[1]/2)
function P_synch(p::Real,η::Real)
    
    if η > 1e-9
        return solve(IntegralProblem(synch_integrand,(0,1),(p,η)),QuadGKJL();abstol=1e-3,reltol=1e-3).u
    else
        return NaN
    end
end


#Janky way of using P(p,η) from Tsunetoe 2025:
const P_synch_const = Float32(inv(P_synch(2.5,0.01)))


function load_textures()
    #No skybox
    return 1,1,1,1
end


"""
Calculate 'global' magnetic field at a position.
    """
@inline function calc_B_field(ray)
    #units are Gauss
    B2 = 0f0
    B3 = 0f0
    B4 = 1f1
    return B2,B3,B4
end


"""
Transform magnetic field from global coordinates to 'local'(but not orthogonal ??) coordinates.
"""
@inline function transform_B_b(B2,B3,B4,g,src_vel)
    b1 = 0.0f0
    
    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,2]*src_vel[i]
    end
    
    b1 += B2*part_sum

    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,3]*src_vel[i]
    end

    b1 += B3*part_sum

    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,4]*src_vel[i]
    end

    b1 += B4*part_sum

    b2 = inv(src_vel[1])*(B2 + b1*src_vel[2])

    b3 = inv(src_vel[1])*(B3 + b1*src_vel[3])

    b4 = inv(src_vel[1])*(B4 + b1*src_vel[4])

    #b^2
    bmag2 = 0.0f0
    for i in 1:4
        bmag2 += g[i,1]*b1^2
    end
    for i in 1:4
        bmag2 += g[i,2]*b2^2
    end
    for i in 1:4
        bmag2 += g[i,3]*b3^2
    end
    for i in 1:4
        bmag2 += g[i,4]*b4^2
    end
    #bmag2 is the square-magnitude of the magnetic field in the fluid frame
    #units of sqrt(bmag2) = Gauss

    return b1,b2,b3,b4,bmag2
end


"""
Calculate angle between magnetic field and ray.
"""
@inline function calc_b_angle(b1,b2,b3,b4,bmag2,ray,g,freq_shift)
    b_dot_k = 0.0f0
    
    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,1]*ray[i+4]
    end
    b_dot_k += b1*part_sum

    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,2]*ray[i+4]
    end
    b_dot_k += b2*part_sum

    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,3]*ray[i+4]
    end
    b_dot_k += b3*part_sum

    part_sum = 0.0f0
    for i in 1:4
        part_sum += g[i,4]*ray[i+4]
    end
    b_dot_k += b4*part_sum

    cos2theta = (b_dot_k)^2/(bmag2*freq_shift^2)

    return cos2theta
end


"""
Calculates various fluid-related parameters for a given position and ray (mutating).
"""
@inline function calc_fluid_params!(fluid_params,Ray,source_vel,g,freq_shift)
    #power law index p
    fluid_params[1] = 2.5f0
    #Fluid frame electron γ_min
    fluid_params[2] = 3f1
    #Fluid frame electron γ_max
    fluid_params[3] = 1f9
    
    B2,B3,B4 = calc_B_field(Ray)
    b1,b2,b3,b4,bmag2 = transform_B_b(B2,B3,B4,g,source_vel)
    #nonthermal electron density
    fluid_params[4] = ne_0*bmag2
    #square cosine of angle between ray and magnetic field
    fluid_params[5] = calc_b_angle(b1,b2,b3,b4,bmag2,Ray,g,freq_shift)
    #cyclotron frequency v_c = 2.8f7 * B Hz
    #v_c = B/(1.07f10 nm) = 9.34f-11 B nm^-1
    #B in gauss
    fluid_params[6] = 9.34f-11 * sqrt(bmag2)
    return nothing
end


function give_n_fluid_params()
    return 6
end


"""
Calculates the spectral radiative coefficients for a BB radiator.
"""
@inline function calc_radiative_coefficients(ray,fluid_params,frequency)
    if is_fire(ray)
        p = fluid_params[1]
        γ_min = fluid_params[2]
        γ_max = fluid_params[3]
        cos2theta = fluid_params[5]
        sintheta = sqrt(1 - cos2theta)

        #base emission scale
        j_ν = ne_0*synch_em*(fluid_params[6]/9.34f-11)
        a_ν = ne_0*synch_ab/frequency

        #electron distribution -> spectrum (Marszewski 2021)
        j_ν *= 3^(p/2)*(p - 1f0)*sintheta/(2*(p + 1f0)*(γ_min^(1f0 - p) - γ_max^(1f0 - p)))
        j_ν *= gamma((3*p - 1f0)/12)*gamma((3*p + 19)/12)*(frequency/(fluid_params[6]*sintheta))^(-(p - 1f0)/2)
        
        a_ν *= 3^((p + 1)/2)*(p - 1f0)/(4*(γ_min^(1f0 - p) - γ_max^(1f0 - p)))
        a_ν *= gamma((3*p + 2)/12)*gamma((3*p + 22)/12)*(frequency/(fluid_params[6]*sintheta))^(-(p + 2f0)/2)
        #anisotropic distribution modification (Tsunetoe 2025)
        Φ = P_synch_const*(1 -0.99f0*cos2theta)^(-p/2)
        j_ν *= Φ
        a_ν *= Φ
        
    else
        j_ν = 0f0
        a_ν = 0f0
    end
    return a_ν,j_ν
end


"""
Determines if there is emmiting material at a location.
"""
@inline function is_fire(position)
    #cylinder of radius 1 centered on z axis
    if sqrt(position[2]^2 + position[3]^2) <= 1
        return true
    else
        return false
    end
end


"""
Keeps the integrator from going too fast.
"""
@inline function pad_max_dt(position,max_dt_scale)
    return max_dt_scale*(1+abs(position[2])+abs(position[3]))
end


"""
Velocity of the material at a position (mutating).
"""
@inline function get_source_velocity!(position,source_vel_buffer)
    source_vel_buffer[1] = 1.0f0
    source_vel_buffer[2] = 0.0f0
    source_vel_buffer[3] = 0.0f0
    source_vel_buffer[4] = 0.0f0
    return nothing
end


"""
Velocity of the material at a position (non-mutating).
"""
function get_source_velocity(position)
    return [1.0f0,0.0f0,0.0f0,0.0f0]
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


function skybox_handling!(ray,raylength::Integer,colors,colors_freq,n_bundle_param,skybox1,skybox1_pix_height,
                          skybox2,skybox2_pix_height)
    #No skybox
    return nothing
end
