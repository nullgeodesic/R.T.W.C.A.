"""
v0.4.1
January 15 2026
Author: Levi Malmström
"""
module BitConsts64

export h,ħ,c,π,k_B,map_scale,no_div_zero,emission_scale,f_factor
export DP_Butcher

#CONSTANTS
#Planck Constant in J/Hz
const h = 6.62607015e-34
#Reduced Planck Constant in J*s
const ħ = h/(2*pi)
#Speed of light in m/s
const c = 299792458
#Boltzman Constant in J/K
const k_B = 1.380649e-23
#How many meters corespond to one unit of the map
const map_scale = Float64(1)
#Protects from dividing by zero in certain situations
const no_div_zero = 1e-24
#2h/c^2, but scaled to ν measured in 1/nm instead of Hz
const emission_scale = 2*h*c*1e27
#h/k_B in nm/K, instead of the normal s/K, because ν is in units of nm^-1, not Hz
const f_factor = h*c*1e9/k_B
#Color match function fit values
#Dormand-Prince Butcher table
const DP_Butcher = Array{Float64}([
[0 0 0 0 0 0 0 0];
[1/5 1/5 0 0 0 0 0 0];
[3/10 3/40 9/40 0 0 0 0 0];
[4/5 44/45 -56/15 32/9 0 0 0 0];
[8/9 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0];
[1 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0];
[1 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]])
end

"""
module BitConsts32

export h,ħ,c,π,k_B,map_scale,no_div_zero,emission_scale,ν_factor
export DP_Butcher

#CONSTANTS
#Planck Constant in J/Hz
const h = Float32(6.62607015e-34)
#Reduced Planck Constant in J*s
const ħ = Float32(h/(2*pi))
#Speed of light in m/s
const c = Int32(299792458)
# π ∈ Irrational doesn't play well with CUDA.jl, so we'll just set it to be type Float32.
const π = Float32(π)
#Boltzman Constant in J/K
const k_B = Float32(1.380649e-23)
#How many meters corespond to one unit of the map
const map_scale = Float32(1)
#Protects from dividing by zero in certain situations
const no_div_zero = Float32(1e-24)
#2h/c^2, but scaled to ν measured in 1/nm instead of Hz
const emission_scale = 2*h*c*1f27
#h/k_B in nm/K, instead of the normal s/K, because ν is in units of nm^-1, not Hz
const ν_factor = h*c*1f9/k_B
#Color match function fit values
#Dormand-Prince Butcher table
const DP_Butcher = Array{Float32}([
[0 0 0 0 0 0 0 0];
[1/5 1/5 0 0 0 0 0 0];
[3/10 3/40 9/40 0 0 0 0 0];
[4/5 44/45 -56/15 32/9 0 0 0 0];
[8/9 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0];
[1 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0];
[1 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
[0 5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]])

end
"""
