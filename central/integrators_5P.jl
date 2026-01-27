"""
v0.4.2
January 27 2026
Author: Levi Malmström
"""


"""
Calculates the derivative of the ray with respect to affine parameter (mutating).
"""
@inline function calc_ray_derivative!(Ray,raylength::Integer,colors_freq,slope,source_vel,g)
    #calculate derivative at a point
    for i in 1:4
        #dx/dl=v
        slope[i]=Ray[i+4]
    end
    for i in 5:8
        #geodesic equation
        a=0
        for j in 5:8
            for k in 5:8
                a -= calc_christoffel_udd(Ray,(i-4,j-4,k-4))*Ray[j]*Ray[k]
            end
        end
        if abs(a) > 1e9
            a=1e9
        end
        slope[i]=a
    end

    #calculate the frequency of the ray in the source frame, by nu=E/hbar, E = -p * u
    get_source_velocity!(Ray,source_vel)
    calc_lower_metric!(Ray,g)
    #using my own four loop instead of Julia's stock matrix multiplication
    #@views freq_shift=-transpose(source_vel)*g*Ray[5:8]
    freq_shift = 0
    for j in 1:4
        part_sum = 0
        for i in 1:4
            part_sum += g[i,j] * Ray[i+4]
        end
        freq_shift -= source_vel[j] * part_sum
    end
    """
    if !isfinite(freq_shift) || freq_shift <= 0
        println("Bad freq shift; ", freq_shift, " ", threadid(), " ",Ray[1:8])
    end
    """
    for i in 9:raylength
        if isodd(i)
            nu=colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the invariant brightness the ray "will" (hence the - sign) gain between here and
            #the camera
            slope[i] = -calc_spectral_emission_coeficient(Ray,nu)*exp(-Ray[i+1])/nu^3
        else
            nu=-colors_freq[ceil(Int,(i-8)/2)]*freq_shift
            #derivative of the optical depth which the ray "will" (hence the - sign) pass through in the "future"
            #(+lambda) to get to the camera
            slope[i]=-calc_spectral_absorbtion_coeficient(Ray,nu)*freq_shift
        end
    end
    return nothing
end


"""
Integrate the ray by one time step using the Runge-Kutta-Dormand-Prince Method.
"""
@inline function RKDP_Step_w_buffer!(Ray,y,last_slope,next_slope,raylength::Integer,
                             stepsize::Real,colors_freq, abs_tol,rel_tol,prev_rejected::Bool, max_dt_scale::Real,
                             buffer,k2,k3,k4,k5,k6,
                             k7,source_vel,g)
    k1 = last_slope
    max_dt = pad_max_dt(Ray, max_dt_scale)

    #check for and deal with singularities
    warn_singularity, temp_stepsize = near_singularity(Ray,stepsize,abs_tol)
    if warn_singularity
        #switch to an euler method scheme that should hopefully dodge the singularity itself
        #Ray -= temp_stepsize*k1
        for i in 1:raylength
            Ray[i] -= temp_stepsize*k1[i]
        end

        keepinbounds!(Ray)
        """
        @views if is_singularity(Ray[1:4]) || find_NaN(Ray[1:8])
            println("Bad Numbers Warning! ", Ray[1:8], " ",threadid())
        end
        """
        for i in 1:8
            if !calc_terminate(Ray,temp_stepsize,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale,2,1)
                calc_ray_derivative!(Ray,raylength,colors_freq,k1,source_vel,g)
                #Ray -= temp_stepsize*k1
                for j in 1:raylength
                    Ray[j] -= temp_stepsize*k1[j]
                end
                keepinbounds!(Ray)
            end
            """
            @views if is_singularity(Ray[1:4])
                println("Singularity Warning! ", Ray[1:8], " ",threadid())
            elseif find_NaN(Ray[1:8])
                println("NaN Warning! ", Ray[1:8], " ",threadid())
            end
            """
        end

        #check for horizons and stuff
        if calc_terminate(Ray,temp_stepsize,colors_freq,raylength,abs_tol,rel_tol,max_dt_scale,2,1)
            #return terminal ray
            return stepsize,false
        end
        #then hand the ray back to the normal Runge-Kutta integrator
        calc_ray_derivative!(Ray,raylength,colors_freq,k1,source_vel,g)

        @views stepsize = -pad_max_dt(Ray, max_dt_scale)
    end

    #'buffer' is used to cut out the memory allocations on array math
    #Calculating k's
    #@views k2 = calc_ray_derivative(Ray .+ k1*stepsize .* DP_Butcher[2,2], raylength,colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*k1[i]*DP_Butcher[2,2]
    end
    #k2 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k2,source_vel,g)

    #@views k3 = calc_ray_derivative(Ray .+ stepsize .* (k1 .* DP_Butcher[3,2] .+ k2 .* DP_Butcher[3,3]),
    #raylength,colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[3,2] + k2[i]*DP_Butcher[3,3])
    end
    #k3 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k3,source_vel,g)

    #@views k4 = calc_ray_derivative(Ray .+ stepsize .* (
    #k1 .* DP_Butcher[4,2] + k2 .* DP_Butcher[4,3] + k3 .*DP_Butcher[4,4]),raylength,colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[4,2] + k2[i]*DP_Butcher[4,3] + k3[i]*DP_Butcher[4,4])
    end
    #k4 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k4,source_vel,g)
    
    #@views k5 = @turbo calc_ray_derivative(Ray .+ stepsize .* (
    #k1 .* DP_Butcher[5,2] .+ k2 .* DP_Butcher[5,3] .+ k3 .* DP_Butcher[5,4] .+ k4 .* DP_Butcher[5,5]),
    #raylength, colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[5,2] + k2[i]*DP_Butcher[5,3] + k3[i]*DP_Butcher[5,4] +
            k4[i]*DP_Butcher[5,5])
    end
    #k5 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k5,source_vel,g)

    #@views k6 = calc_ray_derivative(Ray .+ stepsize .*(k1 .* DP_Butcher[6,2] .+ k2 .* DP_Butcher[6,3] +
    #k3.* DP_Butcher[6,4] .+ k4 .* DP_Butcher[6,5] .+ k5 .* DP_Butcher[6,6]),raylength, colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*(k1[i]*DP_Butcher[6,2] + k2[i]*DP_Butcher[6,3] + k3[i]*DP_Butcher[6,4] +
            k4[i]*DP_Butcher[6,5] + k5[i]*DP_Butcher[6,6])
    end
    #k6 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k6,source_vel,g)

    #@views next_slope=k1 .* DP_Butcher[7,2] .+
    #k3 .* DP_Butcher[7,4] .+ k4 .* DP_Butcher[7,5] .+ k5 .* DP_Butcher[7,6] .+
    #k6 .* DP_Butcher[7,7]
    for i in 1:raylength
        buffer[i] = k1[i]*DP_Butcher[7,2] + k3[i]*DP_Butcher[7,4] + k4[i]*DP_Butcher[7,5] + k5[i]*DP_Butcher[7,6] +
            k6[i]*DP_Butcher[7,7]
    end
    #calculate the derivative to be used in the next calculation
    for i in 1:raylength
        next_slope[i] = buffer[i]
    end
    
    
    #k7 = calc_ray_derivative(Ray .+ stepsize .* next_slope,raylength, colors_freq)
    for i in 1:raylength
        buffer[i] = Ray[i] + stepsize*next_slope[i]
    end
    #k7 = calc_ray_derivative(buffer,raylength,colors_freq)
    calc_ray_derivative!(buffer,raylength,colors_freq,k7,source_vel,g)
    
    #Calculating y's
    #y = Ray .+ stepsize .* next_slope
    
    for i in 1:raylength
        y[i] = Ray[i] + stepsize*next_slope[i]
    end
    

    #Estimate error
    #@views y_hat = Ray .+ stepsize .* (
    #k1 .* DP_Butcher[9,2] .+ k3 .* DP_Butcher[9,4] .+ k4 .* DP_Butcher[9,5] .+ k5 .* DP_Butcher[9,6] .+
    #k6 .* DP_Butcher[9,7] .+ k7 .* DP_Butcher[9,8])
    #delta_y = y-y_hat
    for i in 1:raylength
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
        keepinbounds!(y)
        for i in 1:raylength
            last_slope[i] = next_slope[i]
            Ray[i] = y[i]
        end
        return stepsize,false
    elseif error > 1
        #error too large; send back to manager marked as rejected, with new step size
        if prev_rejected && abs(new_stepsize)>=abs(stepsize)
            #make sure stepsize goes down when the error is too large
            new_stepsize=stepsize/2
        end
        keepinbounds!(Ray)
        return new_stepsize,true
    else
        #error small enough; send updated stats marked as accepted, with new step size
        #keep stepsize from increasing quickly
        if abs(new_stepsize)>5*abs(stepsize)
            new_stepsize=5*stepsize
        end
        keepinbounds!(y)
        for i in 1:raylength
            last_slope[i] = next_slope[i]
            Ray[i] = y[i]
        end
        return new_stepsize,false
    end
end
