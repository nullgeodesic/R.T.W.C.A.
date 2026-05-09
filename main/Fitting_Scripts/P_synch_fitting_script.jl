"""
Author: Levi Malmström
"""

using LsqFit
using Printf

# ── tuneable range ──────────────────────────────────────────────────────────
const P_MIN = 1.0          # must be > 1 (exclusive) in theory; 1.0 is fine numerically
const P_MAX = 5.0
const N_PTS = 500           # number of sample points used for fitting
# ───────────────────────────────────────────────────────────────────────────

#Power Law Synchrotron:
synch_integrand(μ,p) = (1 - 0.99*μ^2)^(-p/2)
function P_synch(p::Real)
    
    return solve(IntegralProblem(synch_integrand,(0,1),p),QuadGKJL();abstol=1e-9,reltol=1e-9).u
end

# ── model function ──────────────────────────────────────────────────────────
#function model(pvec::AbstractVector, params::AbstractVector)
    #C, D, E = params[1], params[2], params[3]
    #return @. C / sqrt(pvec) * 100.0^((pvec - 1) / 2) + D*pvec + E
#end
function model(pvec::AbstractVector, params::AbstractVector)
    C, D, E = params[1], params[2], params[3]
    return @. exp(C + E*pvec)*pvec^(D*pvec)
end

# ── build sample grid and compute reference values ─────────────────────────
p_samples = range(P_MIN, P_MAX; length=N_PTS) |> collect
I_ref     = P_synch.(p_samples)

# ── initial guess for [C, D] ───────────────────────────────────────────────
C0 = 1.0
D0 = 1.0
E0 = 1.0
p0 = [C0, D0, E0]

# ── least-squares fit ──────────────────────────────────────────────────────
fit = curve_fit(model, p_samples, I_ref, p0)
C_fit, D_fit, E_fit = coef(fit)

# ── results ────────────────────────────────────────────────────────────────
println("=" ^ 60)
@printf "Fitted parameters (p ∈ [%.2f, %.2f], %d points):\n" P_MIN P_MAX N_PTS
@printf "  C = %+.6e\n" C_fit
@printf "  D = %+.6e\n" D_fit
@printf "  E = %+.6e\n" E_fit
println()

println("Pointwise comparison (every 5th sample point):")
println("-" ^ 60)
@printf "  %-8s  %-14s  %-14s  %-10s\n" "p" "I(p) exact" "f(p) approx" "rel. error"
println("-" ^ 60)
for i in 1:5:N_PTS
    p_i   = p_samples[i]
    exact = I_ref[i]
    approx = only(model([p_i], [C_fit, D_fit, E_fit]))
    relerr = (approx - exact) / exact * 100
    @printf "  %-8.4f  %-14.6e  %-14.6e  %+.3f%%\n" p_i exact approx relerr
end
println("=" ^ 60)

# ── convenience: evaluate the fitted approximation at any p ────────────────
f_approx(p::Real) = exp(C_fit + E_fit*p)*p^(D_fit*p)
