"""
Author: Levi Malmström
"""

#Dneg Wormhole Hole: (t,l,θ,ϕ);(γ,v,ω,ψ)
#I'm using the parameters from Interstellar's wormhole for this wormhole.
const W_over_M = −log(sec(π/(2*sqrt(2)))) + (π/(2*sqrt(2)))tan(π/(2*sqrt(2)))
#Wormhole radius
const ρ = 1.0f0
#Wormhole half-length
const a_wh = 0.005*ρ
#Lensing Width
const W = 0.05*ρ
const M = W/W_over_M


@inline function x_of_l(l::Real)
    return 2*(abs(l) - a_wh)/(π*M)
end


@inline function r_of_l(l::Real)
    if abs(l) > a_wh
        x = x_of_l(l)
        return ρ + M*(x*atan(x) - log(1+x^2)/2)
    else
        return ρ
    end
end


@inline function dr_dl(l::Real)
    if abs(l) <= a_wh
        return 0.0
    else
        x = x_of_l(l)
        return 2*sign(l)*(atan(x) + (x-1)/(1+x^2))/π
    end
end


"""
Calculates the metric-matrix g_ij at a position (mutating).
"""
@inline function calc_lower_metric!(position,g)
    g[1,1] = -1
    g[2,2] = 1
    g[3,3] = r_of_l(position[2])^2
    g[4,4] = (r_of_l(position[2])*sin(position[3]))^2
    return nothing
end


"""
Calculates the metric-matrix g_ij at a position (non-mutating).
"""
function calc_lower_metric(position)
    g=Matrix{Float64}(I,4,4)
    g[1,1] = -1
    g[2,2] = 1
    g[3,3] = r_of_l(position[2])^2
    g[4,4] = (r_of_l(position[2])*sin(position[3]))^2
    return g
end


"""
Calculates the vierbein at a position.
"""
function calc_vierbein(position)
    vierbein = Matrix{Float64}(I,4,4)
    vierbein[1,1] = 1
    vierbein[2,2] = 1
    vierbein[3,3] = r_of_l(position[2])
    vierbein[4,4] = r_of_l(position[2])*sin(position[3])
    return vierbein
end


"""
Calculates the inverse vierbein at a position.
"""
function calc_inv_vierbein(position)
    inv_vierbein = Matrix{Float64}(I,4,4)
    inv_vierbein[1,1] = 1
    inv_vierbein[2,2] = 1
    inv_vierbein[3,3] = inv(r_of_l(position[2]))
    inv_vierbein[4,4] = inv(r_of_l(position[2])*sin(position[3]) + no_div_zero)
    return inv_vierbein
end


"""
Calculates the Christoffel Symbols of the second kind at a position.
"""
@inline function calc_christoffel_udd(ray,index::Tuple)
    if index[1]==1
        return 0.0
    elseif index[1]==2
        if index[2]==3 && index[3]==3
            return -r_of_l(ray[2])*dr_dl(ray[2])
        elseif index[2]==4 && index[3]==4
            return -r_of_l(ray[2])*dr_dl(ray[2])*sin(ray[3])^2
        else
            return 0.0
        end
    elseif index[1]==3
        if (index[2]==2 && index[3]==3) || (index[2]==3 && index[3]==2)
            return inv(r_of_l(ray[2]))*dr_dl(ray[2])
        elseif index[2]==4 && index[3]==4
            return -sin(ray[3])*cos(ray[3])
        else
            return 0.0
        end
    elseif index[1]==4
        if (index[2]==2 && index[3]==4) || (index[2]==4 && index[3]==2)
            return inv(r_of_l(ray[2]))*dr_dl(ray[2])
        elseif (index[2]==3 && index[3]==4) || (index[2]==4 && index[3]==3)
            return cot(ray[3])
        else
            return 0.0
        end
    else
        return 0.0
    end           
end


"""
Determines if the ray is near a coordinate singularity.
"""
@inline function near_singularity(ray,stepsize::Real,abs_tol)
    dθ = stepsize*ray[7]
    θ_new = ray[3] + dθ

    #check if near θ = 0
    if abs(θ_new) <= abs_tol[3]
        return true, 2*abs(ray[3])/(abs(ray[7]*9) + no_div_zero)
    #check if near θ = π
    elseif abs(θ_new - π) <= abs_tol[3]
        return true, 2*abs((ray[3] - pi)/(abs(ray[7]*9) + no_div_zero))
    else
        return false, stepsize
    end
end


"""
Keeps the ray in the world bounds.
"""
@inline function keepinbounds!(ray)
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
    if abs(position[3]) <= 1e-323 || abs(position[3] - π) <= 1e-323
        return true
    else
        return false
    end
end
