module utils

export C12, C22, Y_eq, z_eq, Mstar_s

using QuadGK
import SpecialFunctions: besselk, besselkx
import PhysicalConstants.CODATA2018: G
using NaturallyUnitful

include((@__DIR__)*"/"*"hubble.jl")
using .Hubble

include((@__DIR__)*"/"*"gstar.jl")
using .gstar

# Physical Constants:
const Grav = ustrip(uconvert(u"GeV^-2", natural(float(G)))) # Gravitational constant
const T0 = 2.348e-13 # GeV, present day temperature
const s0 = ustrip(uconvert(u"GeV^3", natural(2891.2*u"cm^-3"))) # GeV^3, 
                                                # present day entropy density
const Mstar_no_s = 3*(2.131e-42)^2/(2*T0^3*8*π*Grav) # Collection of constants
                                        # for calculating the relic abundance
const lifetime_conv = ustrip(unnatural(u"s", 1*u"GeV^-1")) # conversion from
                                            # GeV^-1 to s, used for lifetime
const h_consts = sqrt(4*π^3*Grav/45) # the constants in the hubble constant
s_factor(T) = (T0^3*2*π^2*gstar_interp(T)/45)/s0 # factor that accounts for the
                                # difference in using Y = n / T^3 and Y = n / s
Mstar_s(T) = Mstar_no_s * s_factor(T)

function bkx_bigx(k,x)
    u = sqrt(1/x)
    t1 = u
    t2 = .5 * (k-.5) * (k+.5) * u^3
    t3 = (1/128) * (2k-3) * (2k-1) * (2k+1) * (2k+3) * u^5
    t4 = (1/3072) * (2k-5) * (2k-3) * (2k-1) * (2k+1) * (2k+3) * (2k+5) * u^7
    return sqrt(π/2) * (t1+t2+t3+t4)
end

function bkx(k,x)
    if x < 1e8
        return besselkx(k,x)
    else
        return bkx_bigx(k,x)
    end
end

"""
    Y_eq(x, g)

Gives the equilibrium yield Y_eq = n_eq / T^3 where n_eq is the equilibrium
number density for a particle X with g spin degrees of freedom.

x : m_X / T, mass of X divided by temperature
g : # of spin degrees of freedom of X
"""
Y_eq(x, g) = g/(2*π^2) * x^2 * besselk(2,x)

"-log(Y_eq)"
z_eq(x,g) = x - log(g*x^2*bkx(2,x)/(2*π^2))

"K_2(x)/K_2(δ*x) expanded at small and large x for numerical stability"
function bk2_ratio(x, δ)
    if x < 1e-2
        return δ^2*(1-(1-δ^2)*x^2/4)
    elseif x < 1e2/δ
        return besselk(2,x)/besselk(2,δ*x)
    else
        return exp(-(1-δ)*x)*(sqrt(δ)-15*(1-δ)/(8*x*sqrt(δ)))
    end
end

"Integral for 2 -> 2 collisions"
function I_22(pi2σ, x; args=())
    quadgk(y -> (y+2*x)^2 * pi2σ(((y+2*x)/x)^2, args...) * bkx(1, y+2*x) * 
            exp(-y), 0, Inf)[1]
end


"""
    C22(σ, x, g, m, args...)

The 2 to 2 collision term that enters into the Boltzmann equation.

σ : dimensionless cross-section (σ = xsec * m^2)
x : m/T
m : mass of particle p that appears in x
g : spin degrees of freedom of particle p
"""
function C22(pi2σ::Function, x, m, g; args=())
    redh = reduced_H(m/x)
    pref = 1 / (x^4*m*redh*g^2*bkx(2,x)^2)
    return pref*I_22(pi2σ, x; args=args)
end

"""
    C12(Γ, x, model)

The 1 to 2 collision term that enters into the Boltzmann equation.

Γ : dimension*full* decay rate [GeV]
x : m/T
m : mass of particle p that appears in x
"""
function C12(Γ::Function, x, m; args=())
    redh = reduced_H(m/x)
    return x*Γ(args...)*bkx(1,x)/(m^2*redh*bkx(2,x))
end

end