"""
Author: Levi Malmström
"""

"""
Calculates the metric-matrix g_ij at a position (mutating).
"""
@inline function calc_lower_metric!(position,g)
    #In this build I'm going with cartesian minkowski space
    g[1,1]=-1
    return nothing
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
Calculates the vierbein at a position.
"""
function calc_vierbein(position)
    #In this build I'm going with cartesian minkowski space
    return Matrix{Float64}(I,4,4)
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
@inline function calc_christoffel_udd(position,index::Tuple)
    #In this build I'm going with cartesian minkowski space
    return 0.0
end


"""
Determines if the ray is near a coordinate singularity.
"""
@inline function near_singularity(ray,stepsize::Real,abs_tol)
    return false, stepsize
end


"""
Keeps the ray in the world bounds.
"""
@inline function keepinbounds!(ray)
    return nothing
end


"""
If the ray is on a coordinate singularity.
"""
@inline function is_singularity(position)
    return false
end
