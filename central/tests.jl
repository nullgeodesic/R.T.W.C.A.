"""
v0.3.0
September 19 2025
Author: Levi Malmström
"""

include("control.jl")

function find_NaN(array)
    has_nan=false
    for n in array
        if isnan(n)
            has_nan=true
        end
    end
    return has_nan
end

function test_rk4(init_ray,raylength,stepsize,colors_freq,n_steps)
    ray=init_ray
    x=zeros(n_steps+1)
    y=zeros(n_steps+1)
    z=zeros(n_steps+1)
    x[1]=init_ray[2]
    y[1]=init_ray[3]
    z[1]=init_ray[4]
    
    for i in 1:n_steps
        #stepsize is negative becaues we're going back in time
        ray=RK4_Step(ray,raylength,-stepsize,colors_freq)        
        x[i+1]=ray[2]
        y[i+1]=ray[3]
        z[i+1]=ray[4]
    end
    return ray,x,y,z
end

function test_cie_functs()
    wavelengths=zeros(Int,351)
    x=zeros(351)
    y=zeros(351)
    z=zeros(351)
    for i in 1:351
        wavelengths[i]=i+374
        x[i]=cie_x(wavelengths[i])
        y[i]=cie_y(wavelengths[i])
        z[i]=cie_z(wavelengths[i])
    end
    plot(wavelengths,x)
    plot!(wavelengths,y)
    plot!(wavelengths,z)
end

function test_planck(min_wavelength,max_wavelength,T)
    colors=range(min_wavelength,step=1,stop=max_wavelength)
    f(x)=c*1e9/x
    colors_freq=f.(colors)

    
    B_nu=calc_planck.(T,colors_freq)
    B_lambda=zeros(length(B_nu))
    for i in 1:length(B_nu)
        B_lambda[i]=B_nu[i]*c/colors[i]^2
    end
    
    plot(colors,B_lambda)
end

function check_spectrum(ray,colors)
    T=get_temp(1)
    fine_colors=range(minimum(colors),step=1,stop=maximum(colors))
    f(x)=c*1e9/x
    colors_freq=f.(colors)
    fine_colors_freq=f.(fine_colors)
    
    I_nu_planck=calc_planck.(T,fine_colors_freq)
    I_planck=zeros(length(I_nu_planck))
    
    I_ray=ray_to_I_lambda(ray,colors,colors_freq)

    for i in 1:length(I_nu_planck)
        I_planck[i]=I_nu_planck[i]*c/fine_colors[i]^2
    end
    I_planck=(maximum(I_ray)/maximum(I_planck))*I_planck
    
    plot(colors,I_ray)
    plot!(fine_colors,I_planck)
end
