"""
Author: Levi Malmström
"""
#CONSTANTS
#Kerr Black Hole in Boyer-Lindquist coordinates: (t,r,θ,ϕ);(γ,v,ω,ψ)
#Mass (c = G = 1) in units of map_scale
const M = 1.0f0
#Schwarzschild Radius
const r_s = 2*M
#r > r_s
#Spin (angular momentum of BH per unit mass)
#const α = 0.998f0 * r_s/2
const α = 0.998f0 * r_s/2
#const α = 0.0f0 * r_s/2
# 0 <= |a| < 1
#event horizon (different for a spinning BH)
const r_ev = r_s/2 + sqrt((r_s^2 / 4) - α^2)


function calc_orb_limits()
    #Spin parameter χ
    χ = 2*α/r_s
    Z_1 = 1 + (1 - χ^2)^(1/3)*((1 + χ)^(3/2) + (1 - χ)^(3/2))
    Z_2 = sqrt(3*χ^2 + Z_1^2)
    #prograde ISCO
    r_ISCO_prograde = Float32(M*(3 + Z_2 - sqrt((3-Z_1)*(3 + Z_1 + 2*Z_2))))
    #retrograde ISCO
    r_ISCO_retrograde = Float32(M*(3 + Z_2 + sqrt((3-Z_1)*(3 + Z_1 + 2*Z_2))))
    #prograde photon orbit
    r_ph_prograde = Float32(r_s*(1 + cos(2*acos(-abs(α)/M)/3)))
    #retrograde photon orbit
    r_ph_retrograde = Float32(r_s*(1 + cos(2*acos(abs(α)/M)/3)))

    return r_ISCO_prograde,r_ISCO_retrograde,r_ph_prograde,r_ph_retrograde
end

const r_ISCO_prograde,r_ISCO_retrograde,r_ph_prograde,r_ph_retrograde = calc_orb_limits()


"""
Calculates values that show up a lot.
"""
@inline function calc_metric_params(position)
    Σ = position[2]^2 + (α*cos(position[3]))^2
    Δ = position[2]^2 - position[2]*r_s + α^2
    A = (position[2]^2 + α^2)*Σ + r_s*position[2]*(α*sin(position[3]))^2
    return Σ,Δ,A
end


"""
Calculates the metric-matrix g_ij at a position (mutating).
"""
@inline function calc_lower_metric!(position,g)
    Σ,Δ,A = calc_metric_params(position)
    g[1,1] = -(1.0f0 - r_s*position[2]/Σ)
    g[1,4] = -2*r_s*α*position[2]*sin(position[3])^2/Σ
    g[4,1] = g[1,4]
    g[2,2] = Σ/Δ
    g[3,3] = Σ
    g[4,4] = (A/Σ)*sin(position[3])^2
    return nothing
end


"""
Calculates the metric-matrix g_ij at a position (non-mutating).
"""
function calc_lower_metric(position)
    g = Matrix{Float64}(I,4,4)
    Σ,Δ,A = calc_metric_params(position)
    g[1,1] = -(1 - r_s*position[2]/Σ)
    g[1,4] = -2*r_s*α*position[2]*sin(position[3])^2/Σ
    g[4,1] = g[1,4]
    g[2,2] = Σ/Δ
    g[3,3] = Σ
    g[4,4] = (A/Σ)*sin(position[3])^2
    return g
end


"""
Calculates the vierbein at a position
"""
function calc_vierbein(position)
    vierbein = Matrix{Float64}(I,4,4)
    Σ,Δ,A = calc_metric_params(position)
    vierbein[1,1] = sqrt(Σ*Δ/A)
    vierbein[4,1] = -(A/r_s*α*position[2])*vierbein[1,1]
    vierbein[2,2] = sqrt(Σ/Δ)
    vierbein[3,3] = sqrt(Σ)
    vierbein[4,4] = sqrt(A/Σ)*sin(position[3])
    return vierbein
end


"""
Calculates the inverse vierbein at a position.
"""
function calc_inv_vierbein(position)
    inv_vierbein = Matrix{Float64}(I,4,4)
    Σ,Δ,A = calc_metric_params(position)
    inv_vierbein[1,1] = sqrt(A/(Σ*Δ))
    inv_vierbein[1,4] = (r_s*α*position[2]/A)*inv_vierbein[1,1]
    inv_vierbein[2,2] = sqrt(Δ/Σ)
    inv_vierbein[3,3] = inv(sqrt(Σ))
    inv_vierbein[4,4] = sqrt(Σ/A)/sin(position[3])
    return inv_vierbein
end


"""
Calculates the Christoffel symbols of the second kind at a position.
"""
@inline function calc_christoffel_udd(ray,index::Tuple)
    t = Float32(ray[1])
    r = Float32(ray[2])
    θ = Float32(ray[3])
    ϕ = Float32(ray[4])
    Σ = r^2 + (α*cos(θ))^2
    Δ = r^2 - r*r_s + α^2
    A = (r^2 + α^2)*Σ + r_s*r*(α*sin(θ))^2
    
    Γ = 0.0f0
    if index[1] == 1
        if (index[2] == 1 && index[3] == 2) || (index[2] == 2 && index[3] == 1)
            Γ = r_s*(r^2 + α^2)*(r^2 - (α*cos(θ))^2)/(2*Σ^2*Δ)
        elseif (index[2] == 1 && index[3] == 3) || (index[2] == 3 && index[3] == 1)
            Γ = -r_s*α^2*r*sin(2*θ)/(2*Σ^2)
        elseif (index[2] == 2 && index[3] == 4) || (index[2] == 4 && index[3] == 2)
            Γ = r_s*α*sin(θ)^2*((α*cos(θ))^2*(α^2 - r^2) - r^2*(α^2 + 3*r^2))/(2*Σ^2*Δ)
        elseif (index[2] == 3 && index[3] == 4) || (index[2] == 4 && index[3] == 3)
            Γ = r_s*α^3*r*sin(θ)^3*cos(θ)/Σ^2
        end
    elseif index[1] == 2
        if index[2] == 1 && index[3] == 1
            Γ = r_s*Δ*(r^2 - α^2*cos(θ)^2)/(2*Σ^3)
        elseif (index[2] == 1 && index[3] == 4) || (index[2] == 4 && index[3] == 1)
            Γ = -Δ*r_s*α*sin(θ)^2*(r^2 - (α*cos(θ))^2)/(2*Σ^3)
        elseif index[2] == 2 && index[3] == 2
            Γ = (2*r*(α*sin(θ))^2 - r_s*(r^2 - (α*cos(θ))^2))/(2*Σ*Δ)
        elseif (index[2] == 2 && index[3] == 3) || (index[2] == 3 && index[3] == 2)
            Γ = -α^2*sin(2*θ)/(2*Σ)
        elseif index[2] == 3 && index[3] == 3
            Γ = -r*Δ/Σ
        elseif index[2] == 4 && index[3] == 4
            Γ = Δ*sin(θ)^2*(-2*r*Σ^2 + r_s*(α*sin(θ))^2*(r^2 - (α*cos(θ))^2))/(2*Σ^3)
        end
    elseif index[1] == 3
        if index[2] == 1 && index[3] == 1
            Γ = -r_s*α^2*r*sin(2*θ)/(2*Σ^3)
        elseif (index[2] == 1 && index[3] == 4) || (index[2] == 4 && index[3] == 1)
            Γ = r_s*α*r*(r^2 + α^2)*sin(2*θ)/(2*Σ^3)
        elseif index[2] == 2 && index[3] == 2
            Γ = α^2*sin(2*θ)/(2*Σ*Δ)
        elseif (index[2] == 2 && index[3] == 3) || (index[2] == 3 && index[3] == 2)
            Γ = r/Σ
        elseif index[2] == 3 && index[3] == 3
            Γ = -α^2*sin(2*θ)/(2*Σ)
        elseif index[2] == 4 && index[3] == 4
            Γ = -sin(2*θ)*(A*Σ + (r^2 + α^2)*r_s*r*(α*sin(θ))^2)/(2*Σ^3)
        end
    elseif index[1] == 4
        if (index[2] == 1 && index[3] == 2) || (index[2] == 2 && index[3] == 1)
            Γ = r_s*α*(r^2 - (α*cos(θ))^2)/(2*Σ^2*Δ)
        elseif (index[2] == 1 && index[3] == 3) || (index[2] == 3 && index[3] == 1)
            Γ = -r_s*α*r*cot(θ)/Σ^2
        elseif (index[2] == 2 && index[3] == 4) || (index[2] == 4 && index[3] == 2)
            Γ = (2*r*Σ^2 + r_s*((α^2*sin(2*θ)/2)^2 - r^2*(Σ + r^2 + α^2)))/(2*Σ^2*Δ)
        elseif (index[2] == 3 && index[3] == 4) || (index[2] == 4 && index[3] == 3)
            Γ = cot(θ)*(Σ^2 + r_s*r*(α*sin(θ))^2)/Σ^2
        end
    end

    return Γ
end


"""
Determines if the ray is near a coordinate singularity.
"""
@inline function near_singularity(ray,stepsize::Real,abs_tol)
    dθ = stepsize*ray[7]
    θ_new = ray[3] + dθ
    dr = stepsize*ray[6]
    new_r = ray[2] + dr

    #check if near θ = 0
    if abs(θ_new) <= abs_tol[3]
        return true, 2*abs(ray[3])/(abs(ray[7]*9) + no_div_zero)
    #check if near θ = π
    elseif abs(θ_new - π) <= abs_tol[3]
        return true, 2*abs((ray[3] - π)/(abs(ray[7]*9) + no_div_zero))
"""
    #check if near r = 0
    elseif abs(new_r) <= abs_tol[2]
        return true, 2*(abs_tol[2])/(abs(ray[6]*9) + no_div_zero)
"""
    else
        return false, stepsize
    end
end


"""
Keeps the ray in the world bounds.
"""
@inline function keepinbounds!(ray)
"""
    if ray[2] <= 0
        ray[2] = -ray[2] + no_div_zero
        ray[6] = -ray[6]
        ray[3] -= π
        ray[7] = -ray[7]
    end
"""
    if ray[3] < 0
        ray[3] = -ray[3]
        ray[4] -= π
        ray[7] = -ray[7]
    end

    ray[3] = mod(ray[3],2π)
    if ray[3] >= π
        ray[3] = 2π - ray[3]
        ray[4] += π
        ray[7] = -ray[7]
    end

    #uncomment the next line to restrict ϕ to [0,2π]
    #ray[4] = mod(ray[4],2π)
    
    return nothing
end


"""
If the ray is on a coordinate singularity.
"""
@inline function is_singularity(position)
    if abs(position[2]) <= 1f-32
        return true
    elseif abs(position[3]) <= 1f-32 || abs(position[3] - π) <= 1e-32
        return true
    else
        return false
    end
end
