module LLFDM

export LFDM, Γ_φm_to_φ0eν, Γ_φ_to_χl, pi2σ_φpφm_to_γγ, pi2σ_φpφm_to_γZ,
    pi2σ_φ0φ0_to_ZZ, pi2σ_φpφm_to_ZZ, pi2σ_φpφm_to_WW, pi2σ_φmφ0_to_WmZ,
    pi2σ_φmφ0_to_Wmγ, pi2σ_φpφm_to_ll, pi2σ_φ0φ0_to_ll, pi2σ_φmφ0_to_lvl,
    pi2σ_φ0φ0_to_WW, Γ_φm_to_φ0μν

using QuadGK
import PhysicalConstants.CODATA2018: FineStructureConstant

include((@__DIR__) * "/../../muc_venv/hubble.jl")
using .Hubble

const αEM = float(FineStructureConstant)
const sw = sqrt(.23121) # sine of weak-mixing angle squared
const cw = sqrt(1 - sw^2)
const mZ = 91.1876 # GeV
const mW = 80.370 # GeV
const mμ = .1056583755 # GeV
const α2 = αEM/sw^2

function arctanh(x)
    .5*(log(1+x) - log(1-x))
end

function arccoth(x)
    .5*(log(x+1) - log(x-1))
end

function _A_f_ΔM(r)
    .5*(r^2 - 2.0 - r*sqrt(r^2-4.0))
end

function _f_ΔM_sca(r)
    -0.25*r*(2*r^3*log(r)+(r^2-4.0)^(3/2)*log(_A_f_ΔM(r)))
end

function deltaM(M, Q, Qp, Y)
    α2*M/(4*π) * ( (Q^2-Qp^2)*sw^2*_f_ΔM_sca(mZ/M) + (Q-Qp)*(Q+Qp-2*Y)*
        (_f_ΔM_sca(mW/M)-_f_ΔM_sca(mZ/M)) )
end

"""
Struct that includes relevant parameters for the one flavor model.
Most observables take this as an input.

λ : Yukawa coupling of operator φ \bar{l} χ

mφ : Mass of new charged scalar φ

gφ : spin degrees of freedom for φ, should be 1

mχ : Mass of new fermion χ

gχ : spin degrees of freedom for χ, should be 2

h : Hubble evaluated at mφ
"""
struct LFDM{T <: Real}
    λ::T
    mφ::T
    gφ::T
    mχ::T
    gχ::T
    δm::T
    h::T
end

function LFDM(λ::Real, mφ::Real, mχ::Real)
    hub = hubble(mφ)
    dm = real(deltaM(complex(mφ), -1, 0, -0.5))
    LFDM(promote(λ, mφ, 1, mχ, 2, dm, hub)...)
end

function LFDM(λ::Real, mφ::Real, gφ::Real, mχ::Real, gχ::Real)
    hub = hubble(mφ)
    dm = real(deltaM(complex(mφ), -1, 0, -0.5))
    LFDM(promote(λ, mφ, gφ, mχ, gχ, dm, hub)...)
end

function Γ_φm_to_φ0eν(model::LFDM)
    αEM^2*model.δm^5 / (30*π*sw^4*mW^4)
end

function Γ_φm_to_φ0μν(model::LFDM)
    if model.δm <= mμ
        return 0.0
    else
        me = mμ/model.mφ
        y = model.δm / model.mφ
        msca = model.mφ
        function f(u)
            (αEM^2*((msca^4*(-1 + u)*y*(-2 + y + u*y)*
                sqrt(me^4 + u^2*y^2*(-2 + u*y)^2 - 2*me^2*(2 - 2*u*y + u^2*y^2))*
                (-(mW^2*(me^2 - 4*(-1 + u*y)^2)) + 
                    msca^2*(me^4 - 4*me^2*(-1 + u*y)^2 + 
                    4*(-1 + u)*u*y^2*(4 - 2*(1 + 2*u)*y + u*(1 + u)*y^2))))/
                (-(mW^4*(-1 + u*y)^2) - me^2*msca^4*
                (me^2*(-1 + y)^2 - (-1 + u)*(-2 + y)*y^2*(-2 + y + u*y)) + 
                msca^2*mW^2*(-((-1 + u)*u*y^2*(4 - 2*(1 + 2*u)*y + u*(1 + u)*y^2)) + 
                    me^2*(2 - 2*(1 + u)*y + (1 + u^2)*y^2))) - 
            (-(me^2*msca^2) + 4*(msca - msca*u*y)^2)*
                log((-4*mW^2*(msca - msca*u*y)^2 + msca^4*(me^2 + (-2 + y)*y)^2 - 
                    msca^4*(2*(-1 + u)*y - (-1 + u^2)*y^2 + 
                        sqrt(me^4 + u^2*y^2*(-2 + u*y)^2 - 2*me^2*(2 - 2*u*y + u^2*y^2)))^2)/
                (-4*mW^2*(msca - msca*u*y)^2 + msca^4*(me^2 + (-2 + y)*y)^2 - 
                    msca^4*(-2*(-1 + u)*y + (-1 + u^2)*y^2 + 
                        sqrt(me^4 + u^2*y^2*(-2 + u*y)^2 - 2*me^2*(2 - 2*u*y + u^2*y^2)))^2))))/
            (64.0*msca^3*π*sw^4)
        end
        return 2*msca^2*y*quadgk(f,me/y,1,rtol=1e-6)[1]
    end
end

function Γ_φ_to_χl(model::LFDM)
    model.λ^2*model.mφ*(1-model.mχ^2/model.mφ^2)^2 / (16*π)
end

function pi2σ_φpφm_to_γγ(x, model::LFDM; approx=true)
    if x < 1e8 || !(approx)
        return 4*π*αEM^2/x^2 * (sqrt((x-4)*x^3) + 2*sqrt((x-4)*x) + 
            2*(x-2)*log( (sqrt(x) - sqrt(x-4)) / (sqrt(x) + sqrt(x-4)) ))
    else
        return 4π*αEM^2 - 8*π*αEM^2*log(x)/x
    end
end

function pi2σ_φpφm_to_γZ(x, model::LFDM; approx=true)
    z = mZ / model.mφ
    if x < 1e6 || !(approx)
        pref = αEM^2*π*(1-2*sw^2)^2 / (24*cw^2*sw^2*x^(3/2)*(x-z^2))
        numer = 2*sqrt(x-4)*(x^2+x*(4-3*z^2)+z^4) - (x-2)*sqrt(x)*(z^2-4)*log(
                (x - sqrt((x-4)*x) - 2) / (x + sqrt((x-4)*x) - 2) )
        return 6*pref*numer
    else
        pref = αEM^2*π*(1-2*sw^2)^2 / (12*cw^2*sw^2)
        t2 = ((z^2-4)*log(x) - 2*z^2 + 2)/x
        return 6*pref*(1 + t2)
    end
end

function pi2σ_φ0φ0_to_ZZ(x, model::LFDM; approx=true)
    z = mZ / model.mφ
    if x < 1e10 || !(approx)
        kinef = sqrt( (x-4) * (x-4*z^2) )
        pref = αEM^2*π / (288*cw^4*sw^4*x)
        t1 = kinef * (4*x + 5*z^4 - 24*z^2 + 16) / (x + z^4 - 4*z^2)
        t2 = -(4*x*(z^2-2)+(z^2-4)^2) * log( (kinef - x + 2*z^2)^2 / 
                (kinef + x - 2*z^2)^2 ) / (x - 2*z^2)
        return 9*pref*(t1+t2)
    else
        return 9*π*αEM^2 / (72*cw^4*sw^4)
    end
end

function pi2σ_φpφm_to_ZZ(x, model::LFDM; approx=true)
    pi2σ_φ0φ0_to_ZZ(x, model::LFDM; approx=approx)*(1-2*sw^2)^4
end

function pi2σ_φpφm_to_WW(x, model::LFDM; approx=true)
    w = mW / model.mφ

    if x < 1e16 || !(approx)
        res = (π*αEM^2*((sqrt((-4 + x)*(-4*w^2 + x))*
                (3*x^2*(4*w^8 + w^6*(48 - 52*x) + x^2*(116 + 11*x) - 
                8*w^2*x*(32 + 15*x) + w^4*(-256 + 292*x + 15*x^2)) + 
                8*sw^2*x*(12*w^10 - 8*x^3*(11 + x) + 16*w^8*(-11 + 2*x) + 
                2*w^2*x^2*(38 + 43*x) + w^4*x*(400 - 152*x - 11*x^2) + 
                w^6*(512 - 248*x + 33*x^2)) - 
                4*sw^4*(48*w^10*(-4 + x) - 8*x^4*(11 + x) + 
                16*w^2*x^3*(-2 + 5*x) + w^4*x^2*(304 - 44*x - 11*x^2) + 
                8*w^6*x*(208 - 64*x + 3*x^2) + w^8*(768 - 656*x + 92*x^2))))/
                ((-4*w^2 + w^4 + x)*(w^2 - cw^2*x)^2) + 
                (12*x*(x*(-2*w^4 + w^2*(8 - 7*x) + 16*x) + 
                4*sw^2*(2*w^6 + 2*w^2*(-4 + x)*x - 4*x^2 + w^4*(-8 + 3*x)))*
                log((-2*w^2 + x - sqrt((-4 + x)*(-4*w^2 + x)))/
                (-2*w^2 + x + sqrt((-4 + x)*(-4*w^2 + x)))))/(-w^2 + cw^2*x)))/
                (1728*sw^4*x^3)

        return 9*res
    else
        return 9*π*αEM^2 / (1728*cw^4*sw^4) * (32*sw^4 - 64*sw^2 + 33)
    end
end

function pi2σ_φ0φ0_to_WW(x, model::LFDM; approx=true)
    w = mW / model.mφ

    if x < 1e16 || !(approx)
        res = (π*αEM^2*((sqrt((-4+x)*(-4*w^2+x))*
                (12*(1+4*sw^2)*w^8+
                4*w^6*(36+48*sw^4+44*sw^2*(-4+x)-39*x)+
                x^2*(348+33*x+12*sw^4*(28+3*x)-4*sw^2*(172+17*x))-
                8*w^2*x*(96+45*x+6*sw^4*(8+5*x)-sw^2*(136+77*x))+
                w^4*(-768+876*x+45*x^2+48*sw^4*(-16+3*x+x^2)-
                4*sw^2*(-512+264*x+23*x^2))))/(-4*w^2+w^4+x)+
                12*(w^2-cw^2*x)*((-2+4*sw^2)*w^4-16*(-1+sw^2)*x+
                w^2*(8-7*x+2*sw^2*(-8+3*x)))*
                log(-((-2*w^2+x+sqrt((-4+x)*(-4*w^2+x)))/
                (2*w^2-x+sqrt((-4+x)*(-4*w^2+x)))))))/
                (1728*sw^4*x*(w^2-cw^2*x)^2)

        return 9*res
    else
        return 9*π*αEM^2 / (1728*cw^4*sw^4) * (36*sw^4 - 68*sw^2 + 33)
    end
end

function pi2σ_φmφ0_to_WmZ(x, model::LFDM; approx=true)
    w = mW / model.mφ
    z = mZ / model.mφ

    if x < 1e16 || !(approx)
        res = (π*αEM^2*(((-1+sw^2)^2*(8*w^2+z^2)*
                (x^1.5+sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))-sqrt(x)*(2+w^2+z^2))
                ^3)/(w^3-w*x)^2+
                (12*(-4+w^2)*x^2*(-4+z^2))/
                (-(sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2)))+sqrt(x)*(w^2-x+z^2))
                +(12*(-4+w^2)*(x-2*sw^2*x)^2*(-4+z^2))/
                (-(sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2)))+sqrt(x)*(w^2-x+z^2))
                -(12*(-4+w^2)*x^2*(-4+z^2))/
                (sqrt(-w^2+x-2*w*z-z^2)*sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(w^2-x+z^2))-
                (12*(-4+w^2)*(x-2*sw^2*x)^2*(-4+z^2))/
                (sqrt(-w^2+x-2*w*z-z^2)*sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(w^2-x+z^2))+
                ((-1+sw^2)^2*(8*w^2+z^2)*
                (sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(2+w^2-x+z^2))^3)/(w^3-w*x)^2+
                (3*(-1+sw^2)*sqrt(x)*
                (sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))-
                sqrt(x)*(2+w^2-x+z^2))^2*
                (-8*w^4+(-1+sw^2)*z^2*(2-x+z^2)+
                w^2*(-16+8*x-9*z^2+sw^2*(16+9*z^2))))/(w^3-w*x)^2-
                (3*(-1+sw^2)*sqrt(x)*
                (sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(2+w^2-x+z^2))^2*
                (-8*w^4+(-1+sw^2)*z^2*(2-x+z^2)+
                w^2*(-16+8*x-9*z^2+sw^2*(16+9*z^2))))/(w^3-w*x)^2-
                (3*x*(-(sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2)))+
                sqrt(x)*(2+w^2-x+z^2))*
                (8*w^6+(-1+sw^2)^2*z^2*
                (4+x^2-12*z^2+z^4+2*x*(-2+z^2))-
                w^4*(-32+32*x-9*z^2+sw^2*(32-32*x+2*z^2)+
                sw^4*(16*x+7*z^2))+
                2*w^2*(16+4*x^2+42*z^2+5*z^4+x*(16-11*z^2)-
                2*sw^2*(16+42*z^2+5*z^4+x*(24-7*z^2))+
                sw^4*(16+42*z^2+5*z^4+x*(32-3*z^2)))))/(w^3-w*x)^2+
                (3*x*(sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(2+w^2-x+z^2))*
                (8*w^6+(-1+sw^2)^2*z^2*
                (4+x^2-12*z^2+z^4+2*x*(-2+z^2))-
                w^4*(-32+32*x-9*z^2+sw^2*(32-32*x+2*z^2)+
                sw^4*(16*x+7*z^2))+
                2*w^2*(16+4*x^2+42*z^2+5*z^4+x*(16-11*z^2)-
                2*sw^2*(16+42*z^2+5*z^4+x*(24-7*z^2))+
                sw^4*(16+42*z^2+5*z^4+x*(32-3*z^2)))))/(w^3-w*x)^2+
                (12*x^1.5*(-((-1+2*sw^2)*
                ((-1+sw^2)*z^4*(-8+z^2)+
                w^4*(-4+(-2+4*sw^2)*x+sw^2*z^2)+
                2*x*(-4+z^2)*(2+(-1+sw^2)*z^2)+
                x^2*(-8+z^2+sw^2*(16-3*z^2))+
                w^2*((2-4*sw^2)*x^2+
                2*(-4+z^2)*(-2+(-1+sw^2)*z^2)+
                x*(12-5*z^2+2*sw^2*(-8+3*z^2))))*
                log((-(sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2)))+
                sqrt(x)*(w^2-x+z^2))/
                (sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(w^2-x+z^2))))+
                (8*x*(2-4*sw^2+x)-x*(12+sw^2*(-16+x)+x)*z^2-
                2*(-1+sw^2)*(-4+x)*z^4-(-1+sw^2)*z^6+
                w^4*(2*(2+x)+sw^2*(-8+z^2))-
                w^2*(16+2*x^2+4*z^2-2*z^4+x*(12-5*z^2)+
                2*sw^2*(-16+z^4+2*x*(-2+z^2))))*
                log(-((sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2))+
                sqrt(x)*(w^2-x+z^2))/
                (sqrt(x)*(-w^2+x-z^2)+
                sqrt(-w^2+x-2*w*z-z^2)*
                sqrt((-4+x)*(-w^2+x+2*w*z-z^2)))))))/
                ((w^2-x)*(w^2-x+z^2))))/(1728*cw^2*sw^4*x^2.5)

        return 9*res
    else
        numer = 32*sw^4*w^2 + sw^4*z^2 - 16*sw^2*w^2 - 2*sw^2*z^2 + 8*w^2 + z^2
        return 9*π*αEM^2 / (864*cw^2*sw^2*w^2) * numer
    end
end

function pi2σ_φmφ0_to_Wmγ(x, model::LFDM; approx=true)
    w = mW / model.mφ

    if x < 1e16 || !(approx)
        res = (π*αEM^2*(2*sqrt(-4+x)*(w^4*(-1+x)+w^2*(2-5*x)*x+x^2*(11+x))+
                3*(-4+w^2)*x^2.5*log(-((sqrt(-4+x)+sqrt(x))/
                (sqrt(-4+x)-sqrt(x))))))/(36*sw^2*x^2.5*(-w^2+x))

        return 6*res
    else
        return 6*π*αEM^2 / (18*sw^2)
    end
end

function pi2σ_φpφm_to_ll(x, model::LFDM)
    z = mZ / model.mφ

    4*(π*((-4+x)/x)^1.5*((1+4*sw^4)*x^2+8*sw^2*(-1-sw^2+2*sw^4)*x*z^2+
    32*sw^4*(-1+sw^2)^2*z^4)*αEM^2)/(384*cw^4*sw^4*(x-z^2)^2)
end

function pi2σ_φ0φ0_to_ll(x, model::LFDM)
    z = mZ / model.mφ
    
    4*π*αEM^2/(384*cw^4*sw^4*(x-z^2)^2)*(8*sw^4-4*sw^2+1)*(x-4)^(1.5)*sqrt(x)
end

function pi2σ_φmφ0_to_lvl(x, model::LFDM)
    w = mW / model.mφ

    4*π*αEM^2*(x-4)^(1.5)*sqrt(x) / (96*sw^4*(x-w^2)^2)
end

end