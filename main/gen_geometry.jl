"""
Author: Levi Malmström
"""

"""
Generates a rotation matrix to rotate the direction the camera is pointing
relative to it's direction of motion.
"""
function gen_intrinsic_rotation_matrix(pointing=[0,0,0])
    #generate rotation matrix
    #'roll'
    R_x=zeros(4,4)
    R_x[1,1]=1
    R_x[2,2]=1
    R_x[3,3]=cos(pointing[1])
    R_x[4,3]=sin(pointing[1])
    R_x[3,4]=-sin(pointing[1])
    R_x[4,4]=cos(pointing[1])

    #'pitch'
    R_y=zeros(4,4)
    R_y[1,1]=1
    R_y[2,2]=cos(pointing[2])
    R_y[4,2]=-sin(pointing[2])
    R_y[3,3]=1
    R_y[2,4]=sin(pointing[2])
    R_y[4,4]=cos(pointing[2])

    #'yaw'
    R_z=zeros(4,4)
    R_z[1,1]=1
    R_z[2,2]=cos(pointing[3])
    R_z[3,2]=sin(pointing[3])
    R_z[2,3]=-sin(pointing[3])
    R_z[3,3]=cos(pointing[3])
    R_z[4,4]=1

    R=R_z*R_y*R_x
    return R
end


"""
Generates a rotation matrix to rotate the direction the camera is moving
relative to the local coordinate tetrad.
"""
function gen_final_rotation_matrix(direction=[0,0])
    #turn angle coords into cartesion direction
    d=[sin(direction[1])*cos(direction[2]),
                             sin(direction[1])*sin(direction[2]),
       cos(direction[1])]
    
    #perpendicular vector for rotation
    k=cross([0,0,1],d)
    mag_k = sqrt(transpose(k)*k)
    if mag_k > 0
        k = k/mag_k
        
        #cross product matrix
        K = zeros(3,3)
        K[1,2] = -k[3]
        K[1,3] = k[2]
        K[2,1] = k[3]
        K[2,3] = -k[1]
        K[3,1] = -k[2]
        K[3,2] = k[1]
        
        #make a 3d rotation matrix
        R_3d = Matrix{Float64}(I,3,3)
        R_3d = R_3d + sin(direction[1])*K + (1-cos(direction[1]))*K*K
        
        #make it 4d
        R = zeros(4,4)
        R[1,1] = 1
        R[2:4,2:4] = R_3d
        #transpose to give proper handedness
        R = transpose(R)
        #compensate for unintended roll induced by the transformation
        #plus a factor of π/2 for some reason I forgot.
        R = R*gen_intrinsic_rotation_matrix([0,0,π/2 + direction[2]])
    else
        #println("Final Rotation: Camera velocity already aligned with FIDO z-axis; returning identity matrix")
        #plus a factor of π/2 for some reason I forgot.
        R = Matrix{Float64}(I,4,4)*gen_intrinsic_rotation_matrix([0,0,π/2 + direction[2]])
    end
    return R
end


"""
Lorentz boost on z axis by velocity β.
"""
function lorentz_boost_z(beta)
    Lambda=zeros(4,4)
    gamma=1/sqrt(1-beta^2)
    Lambda[1,1]=gamma
    Lambda[4,1]=-beta*gamma
    Lambda[1,4]=-beta*gamma
    Lambda[2,2]=1
    Lambda[3,3]=1
    Lambda[4,4]=gamma
    return Lambda
end


"""
Returns a matrix that transforms a four vector based at 'position' from spherical coordinates
to cartesian coordinates; (t,r,θ,ϕ) -> (t,x,y,z).
"""
function R_sph_to_cart(position)
    R = zeros(4,4)
    R[1,1] = 1
    R[2,2] = sin(position[3])*cos(position[4])
    R[2,3] = cos(position[3])*cos(position[4])
    R[2,4] = -sin(position[4])
    R[3,2] = sin(position[3])*sin(position[4])
    R[3,3] = cos(position[3])*sin(position[4])
    R[3,4] = cos(position[4])
    R[4,2] = cos(position[3])
    R[4,3] = -sin(position[3])
    return R
end


"""
Keeping cubemap coordinates on the cube map.
"""
function cube_wrap(i1::Integer,j1::Integer,F1::Integer,x::Integer,y::Integer,N::Integer)
    i2 = i1
    j2 = j1
    F2 = F1
    if !(x != 0 && y != 0)
        if F1 == 1
            if x != 0
                if i1 + x > N
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 2
                elseif i1 + x < 1
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 4
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 1
                end
            else
                if j1 + y > N
                    i2 = mod1(j1 + y,N)
                    j2 = N + 1 - i1
                    F2 = 6
                elseif j1 + y < 1
                    i2 = 1 - j1 - y
                    j2 = i1
                    F2 = 5
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 1
                end
            end
        elseif F1 == 2
            if x != 0
                if i1 + x > N
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 3
                elseif i1 + x < 1
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 1
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 2
                end
            else
                if j1 + y > N
                    i2 = i1
                    j2 = mod1(j1 + y,N)
                    F2 = 6
                elseif j1 + y < 1
                    i2 = i1
                    j2 = mod1(j1 + y, N)
                    F2 = 5
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 2
                end
            end
        elseif F1 == 3
            if x != 0
                if i1 + x > N
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 4
                elseif i1 + x < 1
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 2
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 3
                end
            else
                if j1 + y > N
                    i2 = N - mod1(j1 + y,N) + 1
                    j2 = i1
                    F2 = 6
                elseif j1 + y < 1
                    i2 = N + j1 + y
                    j2 = N - i1 + 1
                    F2 = 5
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 3
                end
            end
        elseif F1 == 4
            if x != 0
                if i1 + x > N
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 1
                elseif i1 + x < 1
                    i2 = mod1(i1 + x,N)
                    j2 = j1
                    F2 = 3
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 4
                end
            else
                if j1 + y > N
                    i2 = N - i1 + 1
                    j2 = N - mod1(j1 + y,N) + 1
                    F2 = 6
                elseif j1 + y < 1
                    i2 = N - i1 + 1
                    j2 = 1 - j1 - y
                    F2 = 5
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 4
                end
            end
        elseif F1 == 5
            if x != 0
                if i1 + x > N
                    i2 = N - j1 + 1
                    j2 = mod1(i1 + x,N)
                    F2 = 3
                elseif i1 + x < 1
                    i2 = j1
                    j2 = 1 - i1 - x
                    F2 = 1
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 5
                end
            else
                if j1 + y > N
                    i2 = i1
                    j2 = mod1(j1 + y,N)
                    F2 = 2
                elseif j1 + y < 1
                    i2 = N - i1 + 1
                    j2 = 1 - j1 - y
                    F2 = 4
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 5
                end
            end
        elseif F1 == 6
            if x != 0
                if i1 + x > N
                    i2 = j1
                    j2 = N - mod1(i1 + x,N) + 1
                    F2 = 3
                elseif i1 + x < 1
                    i2 = N - j1 + 1
                    j2 = N - i1 - x
                    F2 = 1
                else
                    i2 = i1 + x
                    j2 = j1
                    F2 = 6
                end
            else
                if j1 + y > N
                    i2 = N - i1 + 1
                    j2 = N - mod1(j1 + y,N) + 1
                    F2 = 4
                elseif j1 + y < 1
                    i2 = i1
                    j2 = mod1(j1 + y,N)
                    F2 = 2
                else
                    i2 = i1
                    j2 = j1 + y
                    F2 = 6
                end
            end
        end
    else
        #apply shifts in proper order/orientation
        if F1 == 6
            if i1 + x < 1
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,-y,0,N)
            elseif i1 + x > N
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,y,0,N)
            else
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,0,y,N)
            end
        elseif F1 == 5
            if i1 + x > N
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,-y,0,N)
            elseif i1 + x < 1
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,y,0,N)
            else
                i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
                i2,j2,F2 = cube_wrap(i1,j1,F1,0,y,N)
            end
        else
            i1,j1,F1 = cube_wrap(i1,j1,F1,x,0,N)
            i2,j2,F2 = cube_wrap(i1,j1,F1,0,y,N)
        end
    end
    return i2,j2,F2
end


function cubemap_coord_calc(i1::Integer,j1::Integer,F1::Integer,N::Integer)
    i2 = copy(i1)
    j2 = copy(j1)
    
    if F1 == 1
        j2 += N
    elseif F1 == 2
        i2 += N
        j2 += N
    elseif F1 == 3
        i2 += 2*N
        j2 += N
    elseif F1 == 4
        i2 += 3*N
        j2 += N
    elseif F1 == 5
        i2 += N
    else
        i2 += N
        j2 += 2*N
    end

    return i2,j2
end
