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
#jet parameters
#const Ω_jet = 0.5*Ω_horizon = 0.5*α/(r_s*r_ev)
const Ω_jet = 0.5f0*α/(r_s*r_ev)
const Γ_min = 1.2f0
const Γ_max = 1f1
const r_min_nozzle = 10f0
const r_max_nozzle = 75f0
const nozzle_acc = (Γ_max - Γ_min)/(r_max_nozzle - r_min_nozzle)
#Approximation of P_synch(p,0.01), decent accuracy up to p ~ 8.5:
#relative error less than ~17% for p <~6.5
#Fit values:
const C_fit = +7.951834f-01
const D_fit = +9.832989f-01
const E_fit = -5.752283f-01
@inline function Approx_P_synch(p::Real)
    return exp(C_fit + E_fit*p)*p^(D_fit*p)
end


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
    bmag2 = max(1f-12,bmag2)
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

    cos2theta = min((b_dot_k)^2/(bmag2*freq_shift^2),1f0)

    return cos2theta
end


"""
Calculates various fluid-related parameters for a given position and ray (mutating).
"""
@inline function calc_fluid_params!(fluid_params,Ray,source_vel,g,freq_shift)
    #power law index p
    fluid_params[1] = 2.5f0 #* min((1 + abs(Ray[4])/20),4.4f0)
    #Fluid frame electron γ_min
    fluid_params[2] = 3f1
    #Fluid frame electron γ_max
    fluid_params[3] = 1f9
    
    B2,B3,B4 = calc_B_field(Ray)
    b1,b2,b3,b4,bmag2 = transform_B_b(B2,B3,B4,g,source_vel)
    #nonthermal electron density
    fluid_params[4] = ne_0*bmag2
    #square cosine of angle between ray and magnetic field
    fluid_params[5] = calc_b_angle(b1,b2,b3,b4,bmag2,Ray,g,freq_shift) + no_div_zero
    #cyclotron frequency v_c = 2.8f7 * B Hz
    #v_c = B/(1.07f10 nm) = 9.34f-11 B nm^-1
    #B in gauss
    fluid_params[6] = 9.34f-11 * sqrt(bmag2)
    #BB temp of simple disk
    fluid_params[7] = 5778
    return nothing
end


function give_n_fluid_params()
    return 7
end


"""
Calculates the spectral radiative coefficients for a BB radiator.
"""
@inline function calc_radiative_coefficients(ray,fluid_params,frequency)
    zone = region_calc(ray)
    if zone == 1
        #disk
        #a_ν units are 1/map_scale
        a_ν = 1f-2
        j_ν = a_ν*calc_planck(fluid_params[7],frequency)
    elseif zone == 2
        #jet
        p = fluid_params[1]
        γ_min = fluid_params[2]
        γ_max = fluid_params[3]
        cos2theta = fluid_params[5]
        sintheta = sqrt(max(1f-12,1 - cos2theta))

        #base emission scale
        j_ν = ne_0*synch_em*(fluid_params[6]/9.34f-11)
        a_ν = ne_0*synch_ab/frequency

        #electron distribution -> spectrum (Marszewski 2021)
        j_ν *= 3^(p/2)*(p - 1f0)*sintheta/(2*(p + 1f0)*(γ_min^(1f0 - p) - γ_max^(1f0 - p)))
        #if frequency/(fluid_params[6]*sintheta) < 0
            #println("A ", frequency, " ", ray[3], " ",sintheta, " B")
        #end
        j_ν *= gamma((3*p - 1f0)/12)*gamma((3*p + 19)/12)*(frequency/(fluid_params[6]*sintheta))^(-(p - 1f0)/2)
        
        a_ν *= 3^((p + 1)/2)*(p - 1f0)/(4*(γ_min^(1f0 - p) - γ_max^(1f0 - p)))
        a_ν *= gamma((3*p + 2)/12)*gamma((3*p + 22)/12)*(frequency/(fluid_params[6]*sintheta))^(-(p + 2f0)/2)
        #anisotropic distribution modification (Tsunetoe 2025)
        Φ = inv(Approx_P_synch(p))*(1 -0.99f0*cos2theta)^(-p/2)
        j_ν *= Φ
        a_ν *= Φ
    else
        j_ν = 0f0
        a_ν = 0f0
    end
    return a_ν,j_ν
end


"""
What region is the ray in.
"""
@inline function region_calc(position)
    if r_ISCO_prograde <= position[2] <= 12 && 31π/64 <= position[3] <= 33π/64
        #wedge shaped disk in r_ISCO_prograde <= r <= 12, 31π/64 <= θ <= 33π/64
        return 1
    elseif position[2] >= r_min_nozzle && !(0.087263889f0 <= position[3] <= 3.054236111f0) &&
        (1f-3 <= position[3] <= 3.1406f0)
        #cone-shaped jet; r >= r_min_nozzle, !(5 deg <= θ <= 175 deg); vacuum next to poles to discourage
        #numerical issues
        return 2
    else
        #empty space
        return 3
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
    #min speed near poles
    if ray[3] < π/8 || ray[3] > 7π/8
        ϵ = min(Float32(max_dt_scale),1f-1)
    else
        ϵ = Float32(max_dt_scale)
    end
    max_step_size = min(inv(abs(inv(ϵ/(abs(ray[6]))) +
        no_div_zero) + abs(inv(ϵ*min(ray[3],π-ray[3])/(abs(ray[7]) +
        no_div_zero) + no_div_zero)) +
        abs(inv(ϵ/(abs(ray[8])+no_div_zero) + no_div_zero)) + no_div_zero),0.5f0)
    max_step_size = max(max_step_size,1f-12)
    return max_step_size
end


"""
Velocity of the material at a position (mutating).
"""
@inline function get_source_velocity!(position,source_vel_buffer,g)
    zone = region_calc(position)
    if zone == 1
        #prograde equatorial circular orbit (disk)
        Ω = sqrt(M)/(position[2]^(3/2) + α*sqrt(M))
        D = position[2]^(3/2) - 3*M*sqrt(position[2]) + 2*α*sqrt(M)
        ut = (position[2]^(3/2) + α*sqrt(M))/(position[2]^(3/4)*sqrt(D))
        #u^α = u^t*[1,0,0,Ω]
        source_vel_buffer[1] = ut
        source_vel_buffer[2] = 0.0f0
        source_vel_buffer[3] = 0.0f0
        source_vel_buffer[4] = ut*Ω
        return nothing
    elseif zone == 2
        ω = -g[1,4]/(g[4,4] + no_div_zero)
        Ω = min(0.9f0*csc(position[3])/position[2],Ω_jet)
        Σ,Δ,A = calc_metric_params(position)
        #lapse = sqrt(p_thorne^2*Δ/Σ_thorne)
        inv_lapse = sqrt(A/(Σ*Δ))
        #rotating conical jet
        Γ = Γ_min + min(nozzle_acc*(position[2] - r_min_nozzle),Γ_max)
        ut = Γ*inv_lapse
        ur2 = (Γ^2 - 1 - g[4,4]*(Γ*inv_lapse)^2 * (Ω - ω)^2)/g[2,2]
        #if ur2 <= 0
            #println("A: ", ur2, " ", position[2], " ", position[3], " ", g[2,2], " ", g[4,4], " ", inv_lapse," :B")
        #end
        source_vel_buffer[1] = ut
        source_vel_buffer[2] = sqrt((Γ^2 - 1 - g[4,4]*(Γ*inv_lapse)^2 * (Ω - ω)^2)/g[2,2])
        source_vel_buffer[3] = 0f0
        source_vel_buffer[4] = ut*Ω
        return nothing
    else
        #FIDO
        A = sqrt(position[2]^2 + (α*cos(position[3]))^2)
        B = position[2]^2 + α^2 - 2*M*position[2]
        C = (position[2]^2 + α^2)^2 - B*(α*sin(position[3]))^2
        #u^t = 1/lapse rate
        ut = sqrt(C)/(A*sqrt(B))
        Ω = 2*α*M*position[2]/C
        #u^α = ut*[1,0,0,Ω]
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
function get_source_velocity(position,g)
    zone = region_calc(position)
    if zone == 1
        #prograde equatorial circular orbit (disk)
        Ω = sqrt(M)/(position[2]^(3/2) + α*sqrt(M))
        D = position[2]^(3/2) - 3*M*sqrt(position[2]) + 2*α*sqrt(M)
        ut = (position[2]^(3/2) + α*sqrt(M))/(position[2]^(3/4)*sqrt(D))
        #u^α = u^t*[1,0,0,Ω]
        return [ut,0,0,ut*Ω]
    elseif zone == 2
        Σ,Δ,A = calc_metric_params(position)
        #rotating conical jet
        Γ = Γ_min + min(nozzle_acc*(position[2] - r_min_nozzle),Γ_max)
        #β = sqrt(1 - Γ^-2), β_r = β*e_r
        β_r = sqrt((1 - Γ^(-2))*(Δ/Σ))
        ut = inv(sqrt(-(g[1,1] + g[2,2]*β_r^2 + g[3,3]*Ω_jet^2 + 2*g[1,4]*Ω_jet)))
        #u^α = u^t*[1,β_r,0,Ω]
        return [ut,β_r,0,ut*Ω_jet]
    else
        #FIDO
        A = sqrt(position[2]^2 + (α*cos(position[3]))^2)
        B = position[2]^2 + α^2 - 2*M*position[2]
        C = (position[2]^2 + α^2)^2 - B*(α*sin(position[3]))^2
        #u^t = 1/lapse rate
        ut = sqrt(C)/(A*sqrt(B))
        Ω = 2*α*M*position[2]/C
        #u^α = ut*[1,0,0,Ω]
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
