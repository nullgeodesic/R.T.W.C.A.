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
