"""
v0.3.7
December 19 2025
Author: Levi Malmström
"""


"""
Calculates the metric-matrix g_ij at a position (mutating).
"""
function calc_lower_metric!(position,g)
    #In this build I'm going with cartesian minkowski space
    g[1,1]=-1
    return g
end


"""
Calculates the metric-matrix g_ij at a position (non-mutating).
"""
function calc_lower_metric(position)
    #In this build I'm going with cartesian minkowski space
    g=Matrix{Float64}(I,4,4)
    g[1,1]=-1
    return g
end


"""
Calculates the inverse vierbein at a position.
"""
function calc_inv_vierbein(position)
    #In this build I'm going with cartesian minkowski space
    return Matrix{Float64}(I,4,4)
end


"""
Calculates the Christoffel symbols of the second kind at a position.
"""
function calc_christoffel_udd(position,index)
    #In this build I'm going with cartesian minkowski space
    Christoffel=0
    return Christoffel
end


"""
Determines if the ray is near a coordinate singularity.
"""
function near_singularity(ray,stepsize::Real,abs_tol)
    return false, stepsize
end


"""
Keeps the ray in the world bounds.
"""
function keepinbounds!(ray::Vector)
    return nothing
end


"""
If the ray is on a coordinate singularity.
"""
function is_singularity(position)
    return false
end
