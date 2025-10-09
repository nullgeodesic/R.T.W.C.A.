"""
v0.3.2
October 3 2025
Author: Levi Malmström
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

function gen_final_rotation_matrix(direction=[0,0])
    #turn angle coords into cartesion direction
    d=[sin(direction[1])*cos(direction[2]),
                             sin(direction[1])*sin(direction[2]),
       cos(direction[1])]
    
    #perpendicular vector for rotation
    k=cross([0,0,1],d)
    mag_k = sqrt(transpose(k)*k)
    if mag_k > 0
        k=k/mag_k
        
        #cross product matrix
        K=zeros(3,3)
        K[1,2]=-k[3]
        K[1,3]=k[2]
        K[2,1]=k[3]
        K[2,3]=-k[1]
        K[3,1]=-k[2]
        K[3,2]=k[1]
        
        #make a 3d rotation matrix
        R_3d=Matrix{Float64}(I,3,3)
        R_3d=R_3d + sin(direction[1])*K + (1-cos(direction[1]))*K*K
        
        #make it 4d
        R=zeros(4,4)
        R[1,1]=1
        R[2:4,2:4]=R_3d
        #transpose to give proper handedness
        R=transpose(R)
    else
        println("Final Rotation: Camera velocity already aligned with FIDO z-axis; returning identity matrix")
        R=Matrix{Float64}(I,4,4)
    end
    return R
end


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


function calc_lower_metric(position)
    #In this build I'm going with cartesian minkowski space
    g=Matrix{Float64}(I,4,4)
    g[1,1]=-1
    return g
end


function calc_local_vierbein(position)
    #In this build I'm going with cartesian minkowski space
    return Matrix{Float64}(I,4,4)
end


function calc_christoffel_udd(position,index)
    #In this build I'm going with cartesian minkowski space
    Christoffel=0
    return Christoffel
end
