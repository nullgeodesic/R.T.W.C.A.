const Hummer1987 = Array{Float32}([
[8.986940175 -4.009515855 8.808871266e-1 2.640245111e-2 -4.580645915e-2 -3.568055702e-3 2.827798067e-3 3.365860195e-4];
[-8.006936989e-1 9.466021705e-1 9.043402532e-2 -9.608451450e-2 -1.885629865e-2 1.050313890e-2 2.800889961e-3 -1.078209202e-3];
[-3.781305103e-1 1.102726332e-1 -1.543619180e-2 8.310561114e-3 2.179620525e-2 4.259726289e-3 -4.181588794e-3 -1.770208330e-3];
[1.877213132e-2 -1.004885705e-1 -5.483366378e-2 -4.520154409e-3 8.366530426e-3 3.70027390e-3 6.889320423e-4 9.460313195e-5];
[7.300158392e-2 3.576785497e-3 -4.545307025e-3 -1.017965604e-2 -9.530211924e-3 -3.450186162e-3 1.040482914e-3 1.407073544e-3];
[-1.744671550e-3 2.864013856e-2 1.903394837e-2 7.091074494e-3 -9.668371391e-4 -2.999107465e-3 -1.820642230e-3 -3.874082085e-4];
[-1.707268366e-2 -4.694254776e-3 1.311691517e-3 5.316703136e-3 5.178193095e-3 2.451228935e-3 -2.277321615e-5 -8.182359057e-4];
[2.567331664e-4 -9.155339970e-3 -6.997479192e-3 -3.571518641e-3 -2.096101038e-4 1.553822487e-3 1.509584686e-3 6.212627837e-4];
[4.098322531e-3 1.635218463e-3 -5.918883504e-4 -2.333091048e-3 -2.484138313e-3 -1.359996060e-3 -5.371426147e-5 5.553549563e-4];
[3.837562402e-5 2.938325230e-3 2.393747064e-3 1.328839809e-3 9.135013312e-5 -7.137252303e-4 -7.656848158e-4 -3.504683798e-4];
[-8.491991820e-4 -3.615327726e-4 3.148015257e-4 8.909207650e-4 9.869737522e-4 6.134671184e-4 1.068883394e-4 -2.046080100e-4]
])

const Gaunt_LUT = Array{Float32}(DataFrame(load(File(format"CSV","Gaunt_Factors/gauntff_merged_Z01.dat"),spacedelim=true,commentchar='#',header_exists=false,skiplines_begin=42)))

#bremstrahlung emissivity in erg cm^3 s^-1 Hz^-1 sr^-1
const C_Yar_cgs = 7.089314638783696e-44
#bremstrahlung emissivity in W m^3 Hz^-1 sr^-1
const C_Yar_si = 7.089314638783696e-57

#Thermal Bremstrahlung:
"""
Calculate Gaunt factor using Chebyshev polynomial fit from Hummer 1987 (note: DON'T USE! I think I've implemented this wrong).
"""
function gaunt_therm_NR_Chebyshev(gam2::Real,u::Real)
    #gam2 = Ry/(k_B_eV*T)
    #u = f_factor*f/T
    x_u = (2*log10(u) - 2.5)/5.5f0
    x_g = log10(gam2)/3

    g_ff = 0.0f0

    for j in 1:8
        C = 0.0f0
        for i in 1:11
            if i != 1
                C += Hummer1987[i,j]*Chebyshev_r(i-1,x_g)
            else
                C += Hummer1987[i,j]*Chebyshev_r(i-1,x_g)/2
            end
        end
        if j != 1
            g_ff += C*Chebyshev_r(j-1,x_u)
        else
            g_ff += C*Chebyshev_r(j-1,x_u)/2
        end
    end

    return g_ff
 end


function gaunt_ei_therm_van_Hoof(gam2::Real,u::Real)
    #tests on the CPU indicate this takes roughly twice the time of the calc_planck function.
    #gam2 = Ry/(k_B_eV*T)
    #u = h*f_Hz/k_B*T = f_factor*f/T
    loggam = Float32(log10(gam2))
    logu = Float32(log10(u))
    
    frac_i = Float32((logu + 16)/0.2f0)
    frac_j = Float32((loggam + 6)/0.2f0)
    floor_i = floor(Int32,frac_i)
    ceil_i = ceil(Int32,frac_i + 1f-3)
    floor_j = floor(Int32,frac_j)
    ceil_j = ceil(Int32,frac_j + 1f-3)

    if (1 <= floor_i) && (1 <= floor_j) && (ceil_i <= 146) && (ceil_j <= 81)
        #bilinear interpolation
        Gaunt1 = Float32((ceil_i - frac_i)*Gaunt_LUT[floor_i,floor_j] + (frac_i - floor_i)*Gaunt_LUT[ceil_i,floor_j])
        Gaunt2 = Float32((ceil_i - frac_i)*Gaunt_LUT[floor_i,ceil_j] + (frac_i - floor_i)*Gaunt_LUT[ceil_i,ceil_j])
        Gaunt = Float32((ceil_j - frac_j)*Gaunt1 + (frac_j - floor_j)*Gaunt2)
    else
        Gaunt = 1.2f0
    end

    return Gaunt
end

#Thermal Transrelativistic Synchrotron:

#Thermal Relativistic Synchrotron:
"""

"""
function j_ν_therm_rel_synch(f::Real,O::Real,sintheta::Real)
    aaaaaaaaaaaa

#Kappa Relativistic Synchrotron:
