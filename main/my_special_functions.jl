"""
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


"""
2D normalized Gaussian function.
"""
function D2_Gauss(x::Real,y::Real;x_0=0.0::Real,y_0=0.0::Real,σ_x=1.0::Real,σ_y=1.0::Real,θ=0.0::Real)
    a = cos(θ)^2 /(2*σ_x^2) + sin(θ)^2 /(2*σ_y^2)
    b = -sin(θ)*cos(θ) /(2*σ_x^2) + sin(θ)*cos(θ) /(2*σ_y^2)
    c = sin(θ)^2 /(2*σ_x^2) + cos(θ)^2 /(2*σ_y^2)
    return inv(2*π*σ_x*σ_y)*exp(-(a*(x-x_0)^2 + 2*b*(x-x_0)*(y-y_0) + c*(y-y_0)^2))
end


"""
Figures out the new shape of a gaussian after non-uniform scaling of it's inputs by a and b.
"""
function distort_gauss(σ_x::Real,σ_y::Real,θ::Real,a::Real,b::Real)
    #use this weird trick from a chatbot: represent gaussian as covariance matrix and do
    #some math to get this result
    Σ_xx_new = a^2*(σ_x^2*cos(θ)^2 + σ_y^2*sin(θ)^2)
    Σ_yy_new = b^2*(σ_x^2*sin(θ)^2 + σ_y^2*cos(θ)^2)
    Σ_xy_new = a*b*(σ_x^2 - σ_y^2)*sin(2*θ)/2
    
    σ_x_new = sqrt((Σ_xx_new + Σ_yy_new + sqrt((Σ_xx_new - Σ_yy_new)^2 + 4*(Σ_xy_new)^2))/2)
    σ_y_new = sqrt((Σ_xx_new + Σ_yy_new - sqrt((Σ_xx_new - Σ_yy_new)^2 + 4*(Σ_xy_new)^2))/2)
    θ_new = atan(2*Σ_xy_new,Σ_xx_new - Σ_yy_new)/2

    return σ_x_new,σ_y_new,θ_new
end


"""
Type 1 Chebyshev polynomials. (recursive)
"""
function Chebyshev_r(n::Integer,x::Real)
    if n > 1
        T = Float32(2*x*Chebyshev_r(n - 1,x) - Chebyshev_r(n - 2,x))
    elseif n == 1
        T = Float32(x)
    else
        T = 1.0f0
    end
    return T
end

"""
Type 1 Chebyshev polynomials. (first 11)
"""
@inline function Chebyshev_e(n::Integer,x::Real)
    if n <= 0
        return 1.0f0
    elseif n == 1
        return Float32(x)
    elseif n == 2
        return Float32(2*x^2 - 1)
    elseif n == 3
        return Float32(4*x^3 - 3*x)
    elseif n == 4
        return Float32(8*x^4 - 8*x^2 + 1)
    elseif n == 5
        return Float32(16*x^5 - 20*x^3 + 5*x)
    elseif n == 6
        return Float32(32*x^6 - 48*x^4 + 18*x^2 - 1)
    elseif n == 7
        return Float32(64*x^7 - 112*x^5 + 56*x^3 - 7*x)
    elseif n == 8
        return Float32(128*x^8 - 256*x^6 + 160*x^4 -32*x^2 + 1)
    elseif n == 9
        return Float32(256*x^9 -576*x^7 + 432*x^5 -120*x^3 + 9*x)
    elseif n == 10
        return Float32(512*x^10 -1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1)
    else
        return 1.0f0
    end
end
