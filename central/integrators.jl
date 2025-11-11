"""
v0.3.2
October 30 2025
Author: Levi Malmström
"""
function calc_ray_derivative(Ray,raylength,colors_freq)
    #calculate derivative at a point
    slope = similar(Ray)
    for i in 1:4
        #dx/dl=v
        slope[i]=Ray[i+4]
    end
    for i in 5:8
        #geodesic equation
        a=0
        for j in 5:8
            for k in 5:8
                @views a -= calc_christoffel_udd(Ray[1:4],CartesianIndex(i-4,j-4,k-4))*Ray[j]*Ray[k]
            end
        end
        slope[i]=a
    end

    #calculate the frequency of the ray in the source frame, by nu=E/hbar, E = -p * u
    @views freq_shift=-transpose(get_source_velocity(Ray[1:4]))*calc_lower_metric(Ray[1:4])*Ray[5:8]

    for i in 9:raylength
        if isodd(i)
            nu=colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the invariant brightness the ray "will" (hence the - sign) gain between here and the camera
            @views slope[i]=-calc_spectral_emission_coeficient(Ray[1:8],nu)*exp(-Ray[i+1])/nu^3
        else
            nu=-colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the optical depth which the ray "will" (hence the - sign) pass through in the "future" (+lambda) to get to the camera
            @views slope[i]=-calc_spectral_absorbtion_coeficient(Ray[1:8],nu)*freq_shift
        end
    end
    return slope
end

function RK4_Step(Ray,raylength,stepsize,colors_freq)
    #Runge-Kutta Method -- Classic
    k1 = calc_ray_derivative(Ray,raylength,colors_freq)
    k2 = calc_ray_derivative(Ray + k1*stepsize/2, raylength,colors_freq)
    k3 = calc_ray_derivative(Ray + k2*stepsize/2, raylength,colors_freq)
    k4 = calc_ray_derivative(Ray + k3*stepsize, raylength,colors_freq)
    return Ray + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
end

function RKDP_Step!(Ray,slope_last,raylength::Integer,stepsize::Real,colors_freq,abs_tol,rel_tol,prev_rejected::Bool,
                   max_dt::Real)
    #Runge-Kutta Method -- Dormand-Prince method
    #Calculating k's
    k1 = slope_last
    @views k2 = calc_ray_derivative(Ray .+ k1*stepsize .* DP_Butcher[2,2], raylength,colors_freq)
    @views k3 = calc_ray_derivative(Ray .+ stepsize .* (k1 .* DP_Butcher[3,2] .+ k2 .* DP_Butcher[3,3]),
                                    raylength,colors_freq)
    @views k4 = calc_ray_derivative(Ray .+ stepsize .* (
        k1 .* DP_Butcher[4,2] + k2 .* DP_Butcher[4,3] + k3 .*DP_Butcher[4,4]),
                             raylength,colors_freq)
    @views k5 = @turbo calc_ray_derivative(Ray .+ stepsize .* (
        k1 .* DP_Butcher[5,2] .+ k2 .* DP_Butcher[5,3] .+ k3 .* DP_Butcher[5,4] .+ k4 .* DP_Butcher[5,5]),
                                    raylength, colors_freq)
    @views k6 = calc_ray_derivative(Ray .+ stepsize .*(
        k1 .* DP_Butcher[6,2] .+ k2 .* DP_Butcher[6,3] + k3.* DP_Butcher[6,4] .+
        k4 .* DP_Butcher[6,5] .+ k5 .* DP_Butcher[6,6]),raylength, colors_freq)
    
    @views next_slope=k1 .* DP_Butcher[7,2] .+
        k3 .* DP_Butcher[7,4] .+ k4 .* DP_Butcher[7,5] .+ k5 .* DP_Butcher[7,6] .+
        k6 .* DP_Butcher[7,7]
    
    k7 = calc_ray_derivative(Ray .+ stepsize .* next_slope,raylength, colors_freq)
    
    #Calculating y's
    y = Ray .+ stepsize .* next_slope
    @views y_hat = Ray .+ stepsize .* (
        k1 .* DP_Butcher[9,2] .+ k3 .* DP_Butcher[9,4] .+ k4 .* DP_Butcher[9,5] .+ k5 .* DP_Butcher[9,6] .+
        k6 .* DP_Butcher[9,7] .+ k7 .* DP_Butcher[9,8])

    #Estimate error
    delta_y = y-y_hat
    error=0
    for i in 1:raylength
        tol = abs_tol[i]+max(abs(Ray[i]),abs(y[i]))*rel_tol[i]
        error += (delta_y[i]/tol)^2
    end
    error=sqrt(error/raylength)

    """
    if isnan(stepsize)
        println(true)
    end
    """

    #decide next step
    #calc size of next step
    new_stepsize=0.9*stepsize*inv(error)^(1/5)

    if abs(new_stepsize)>max_dt
        new_stepsize=-max_dt
    end
    
    if error == 0 || !isfinite(new_stepsize)
        #error small enough; send updated stats marked as accepted, with same step size (to keep out NaN's)
        return y,next_slope,stepsize,false
    elseif error > 1
        #error too large; send back to manager marked as rejected, with new step size
        if prev_rejected && abs(new_stepsize)>=abs(stepsize)
            #make sure stepsize goes down when the error is too large
            new_stepsize=stepsize/2
        end
        return Ray,slope_last,new_stepsize,true
    else
        #error small enough; send updated stats marked as accepted, with new step size
        #keep stepsize from increasing quickly
        if abs(new_stepsize)>5*abs(stepsize)
            new_stepsize=5*stepsize
        end
        return y,next_slope,new_stepsize,false
    end
end

function RKDP_Step_w_buffer!(Ray::Vector,slope_last::Vector,raylength::Integer,stepsize::Real,colors_freq::Vector,
                             abs_tol,rel_tol,prev_rejected::Bool, max_dt::Real,buffer::Vector)
    #Runge-Kutta Method -- Dormand-Prince method
    #'buffer' is used to cut out the memory allocations on array math
    #Calculating k's
    k1 = slope_last
    next_slope=similar(k1)
    for i in eachindex(buffer)
        buffer[i] = k1[i]
    end

    #@views k2 = calc_ray_derivative(Ray .+ k1*stepsize .* DP_Butcher[2,2], raylength,colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*k1[i]*DP_Butcher[2,2]
    end
    k2 = calc_ray_derivative(buffer,raylength,colors_freq)

    #@views k3 = calc_ray_derivative(Ray .+ stepsize .* (k1 .* DP_Butcher[3,2] .+ k2 .* DP_Butcher[3,3]),
    #raylength,colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[3,2] + k2[i]*DP_Butcher[3,3])
    end
    k3 = calc_ray_derivative(buffer,raylength,colors_freq)

    #@views k4 = calc_ray_derivative(Ray .+ stepsize .* (
    #k1 .* DP_Butcher[4,2] + k2 .* DP_Butcher[4,3] + k3 .*DP_Butcher[4,4]),raylength,colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[4,2] + k2[i]*DP_Butcher[4,3] + k3[i]*DP_Butcher[4,4])
    end
    k4 = calc_ray_derivative(buffer,raylength,colors_freq)
    
    #@views k5 = @turbo calc_ray_derivative(Ray .+ stepsize .* (
    #k1 .* DP_Butcher[5,2] .+ k2 .* DP_Butcher[5,3] .+ k3 .* DP_Butcher[5,4] .+ k4 .* DP_Butcher[5,5]),
    #raylength, colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[5,2] + k2[i]*DP_Butcher[5,3] + k3[i]*DP_Butcher[5,4] +
            k4[i]*DP_Butcher[5,5])
    end
    k5 = calc_ray_derivative(buffer,raylength,colors_freq)

    #@views k6 = calc_ray_derivative(Ray .+ stepsize .*(k1 .* DP_Butcher[6,2] .+ k2 .* DP_Butcher[6,3] +
    #k3.* DP_Butcher[6,4] .+ k4 .* DP_Butcher[6,5] .+ k5 .* DP_Butcher[6,6]),raylength, colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[6,2] + k2[i]*DP_Butcher[6,3] + k3[i]*DP_Butcher[6,4] +
            k4[i]*DP_Butcher[6,5] + k5[i]*DP_Butcher[6,6])
    end
    k6 = calc_ray_derivative(buffer,raylength,colors_freq)

    #@views next_slope=k1 .* DP_Butcher[7,2] .+
    #k3 .* DP_Butcher[7,4] .+ k4 .* DP_Butcher[7,5] .+ k5 .* DP_Butcher[7,6] .+
    #k6 .* DP_Butcher[7,7]
    for i in eachindex(buffer)
        buffer[i] = k1[i]*DP_Butcher[7,2] + k3[i]*DP_Butcher[7,4] + k4[i]*DP_Butcher[7,5] + k5[i]*DP_Butcher[7,6] +
            k6[i]*DP_Butcher[7,7]
    end
    next_slope=copy(buffer)
    
    #k7 = calc_ray_derivative(Ray .+ stepsize .* next_slope,raylength, colors_freq)
    for i in eachindex(buffer)
        buffer[i] = Ray[i] + stepsize*next_slope[i]
    end
    k7 = calc_ray_derivative(buffer,raylength,colors_freq)
    
    #Calculating y's
    y = Ray .+ stepsize .* next_slope
    #@views y_hat = Ray .+ stepsize .* (
    #k1 .* DP_Butcher[9,2] .+ k3 .* DP_Butcher[9,4] .+ k4 .* DP_Butcher[9,5] .+ k5 .* DP_Butcher[9,6] .+
    #k6 .* DP_Butcher[9,7] .+ k7 .* DP_Butcher[9,8])

    #Estimate error
    #delta_y = y-y_hat
    for i in eachindex(buffer)
        buffer[i] = y[i] - Ray[i] - stepsize*(k1[i]*DP_Butcher[9,2] + k3[i]*DP_Butcher[9,4] +
            k4[i]*DP_Butcher[9,5] + k5[i]*DP_Butcher[9,6] + k6[i]*DP_Butcher[9,7] + k7[i]*DP_Butcher[9,8])
    end
    
    
    error=0
    for i in 1:raylength
        tol = abs_tol[i]+max(abs(Ray[i]),abs(y[i]))*rel_tol[i]
        error += (buffer[i]/tol)^2
    end
    error=sqrt(error/raylength)

    #decide next step
    #calc size of next step
    new_stepsize=0.9*stepsize*inv(error)^(1/5)

    if abs(new_stepsize)>max_dt
        new_stepsize=-max_dt
    end
    
    if error == 0 || !isfinite(new_stepsize)
        #error small enough; send updated stats marked as accepted, with same step size (to keep out NaN's)
        return y,next_slope,stepsize,false,buffer
    elseif error > 1
        #error too large; send back to manager marked as rejected, with new step size
        if prev_rejected && abs(new_stepsize)>=abs(stepsize)
            #make sure stepsize goes down when the error is too large
            new_stepsize=stepsize/2
        end
        return Ray,slope_last,new_stepsize,true,buffer
    else
        #error small enough; send updated stats marked as accepted, with new step size
        #keep stepsize from increasing quickly
        if abs(new_stepsize)>5*abs(stepsize)
            new_stepsize=5*stepsize
        end
        return y,next_slope,new_stepsize,false,buffer
    end
end
