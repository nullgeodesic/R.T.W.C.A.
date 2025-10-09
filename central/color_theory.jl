"""
v0.3.2
October 3 2025
Author: Levi Malmström
"""

function Selector_funct(x,y,z)
    #Selector function for the fits of the color matching functions
    return y*(1-Heaviside(x))+z*Heaviside(x)
end

function cie_x(wavelength)
    #A fit for the CIE X color matching function
    x=0
    for i in 1:3
        x=x+cie_matrix[1,i]*exp((-1/2)*((wavelength-cie_matrix[2,i])*Selector_funct(
            wavelength-cie_matrix[2,i],cie_matrix[3,i],cie_matrix[4,i]))^2)
    end
    return x
end

function cie_y(wavelength)
    #A fit for the CIE Y color matching function
    y=0
    for i in 1:2
        y=y+cie_matrix[1,i+3]*exp((-1/2)*((wavelength-cie_matrix[2,i+3])*Selector_funct(
            wavelength-cie_matrix[2,i+3],cie_matrix[3,i+3],cie_matrix[4,i+3]))^2)
    end
    return y
end

function cie_z(wavelength)
    #A fit for the CIE Z color matching function
    z=0
    for i in 1:2
        z=z+cie_matrix[1,i+5]*exp((-1/2)*((wavelength-cie_matrix[2,i+5])*Selector_funct(
            wavelength-cie_matrix[2,i+5],cie_matrix[3,i+5],cie_matrix[4,i+5]))^2)
    end
    return z
end

function ray_to_I_λ(ray,colors,colors_freq)
    #Calculates the I_λs from a ray
    n_colors=length(colors)
    I_λs=zeros(n_colors)
    for i in 1:n_colors
        #I_λ=I_nu *c/λ^2, where c=λ*nu
        #since colors are in nm, this gives the spectral radiance in W sr^-1 m^-2 nm^-1
        I_λs[i]=ray[7+2*i]*colors_freq[i]^3*c/colors[i]^2
    end
    return I_λs
end

function calc_xyY(ray,colors,colors_freq)
    #Calculates the color of a ray
    I_λs=ray_to_I_λ(ray,colors,colors_freq)
    
    I_interpolation=linear_interpolation(colors,I_λs)
    
    delXYZ(λ,p) = I_interpolation(λ)*[cie_x(λ),cie_y(λ),cie_z(λ)]
    domain = (minimum(colors),maximum(colors))
    prob=IntegralProblem(delXYZ,domain)
    CIEXYZ=solve(prob,HCubatureJL();reltol=1e-3,abstol=1e-3)   
    if CIEXYZ[1] != 0 && CIEXYZ[2] != 0 && CIEXYZ[3] != 0
        CIEXYZ=XYZ{Float64}(CIEXYZ[1],CIEXYZ[2],CIEXYZ[3])
        CIExyY=xyY(CIEXYZ)
    else
        CIExyY=xyY{Float64}(1,1,0)
    end
    
    return CIExyY
end
