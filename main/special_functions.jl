"""
v0.4.0
December 24 2025
Author: Levi Malmström
"""


"""
Heaviside step function.
"""
function Heaviside(x)
    if x >= 0
        return 1.0
    else
        return 0.0
    end
end

"""
Normalized Gaussian function.
"""
function Normal_Gauss(x::Real,σ::Real,μ::Real)
    return exp(-(x-μ)^2/(2*σ^2))/(σ*sqrt(2*π))
end
