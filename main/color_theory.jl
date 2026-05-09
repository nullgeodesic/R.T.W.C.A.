"""
Author: Levi Malmström
"""
#Rgb spectral scales
#R = 
#G = 
#B = 


#Color match function fit values
const cie_matrix=[
[0.362 1.056 -0.065 0.821 0.286 1.217 0.681];
[442.0 599.8 501.1 568.8 530.9 437.0 459.0];
[0.0624 0.0264 0.0490 0.0213 0.0613 0.0845 0.0385];
[0.0374 0.0323 0.0382 0.0247 0.0322 0.0278 0.0725]]

"""
Selector function for the fits of the color matching functions.
"""
function Selector_funct(x,y,z)
    return y*(1-Heaviside(x))+z*Heaviside(x)
end


"""
A fit for the CIE X color matching function (~red).
"""
function cie_x(wavelength)
    x=0
    for i in 1:3
        x=x+cie_matrix[1,i]*exp((-1/2)*((wavelength-cie_matrix[2,i])*Selector_funct(
            wavelength-cie_matrix[2,i],cie_matrix[3,i],cie_matrix[4,i]))^2)
    end
    return x
end


"""
A fit for the CIE Y color matching function (~green).
"""
function cie_y(wavelength)
    y=0
    for i in 1:2
        y=y+cie_matrix[1,i+3]*exp((-1/2)*((wavelength-cie_matrix[2,i+3])*Selector_funct(
            wavelength-cie_matrix[2,i+3],cie_matrix[3,i+3],cie_matrix[4,i+3]))^2)
    end
    return y
end


"""
A fit for the CIE Z color matching function (~blue).
"""
function cie_z(wavelength)
    z=0
    for i in 1:2
        z=z+cie_matrix[1,i+5]*exp((-1/2)*((wavelength-cie_matrix[2,i+5])*Selector_funct(
            wavelength-cie_matrix[2,i+5],cie_matrix[3,i+5],cie_matrix[4,i+5]))^2)
    end
    return z
end


"""
Calculates the I_λs from a ray.
"""
function ray_to_I_λ(ray,colors,colors_freq)
    n_colors = length(colors)
    I_λs = zeros(n_colors)
    for i in 1:n_colors
        #I_λ=I_nu *c/λ^2, where c=λ*nu
        #since colors are in nm, this gives the spectral radiance in W sr^-1 m^-2 nm^-1
        I_λs[i] = ray[7+2*i]*colors_freq[i]^3*c/colors[i]^2
    end
    return I_λs
end


function precompute_CIE_weights(colors)
    step = colors[2] - colors[1]
    n = length(colors)
    # Simpson's rule weights: 1/3, 4/3, 2/3, 4/3, ..., 1/3 (times step/2)
    # For simplicity, trapezoidal weights are usually sufficient at 1 nm:
    w = fill(1.0*step, n)
    w[1] /= 2; w[end] /= 2  # trapezoidal endpoint correction

    W_X = [cie_x(colors[i]) * w[i] for i in 1:n]
    W_Y = [cie_y(colors[i]) * w[i] for i in 1:n]
    W_Z = [cie_z(colors[i]) * w[i] for i in 1:n]
    return W_X, W_Y, W_Z
end

"""
Fast and accurate (I think) method
"""
function calc_xyY(ray, colors, colors_freq,W_X,W_Y,W_Z)
    I_λs = ray_to_I_λ(ray, colors, colors_freq)
    
    CIEX = sum(I_λs[i] * W_X[i] for i in 1:length(colors))
    CIEY = sum(I_λs[i] * W_Y[i] for i in 1:length(colors))
    CIEZ = sum(I_λs[i] * W_Z[i] for i in 1:length(colors))
    
    if CIEX != 0 && CIEY != 0 && CIEZ != 0
        return xyY(XYZ{Float64}(CIEX, CIEY, CIEZ))
    else
        return xyY{Float64}(1, 1, 0)
    end
end
