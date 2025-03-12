module QLFlavorDM

export QLFDM, Γ_φu_to_χq, Γ_φu_to_χt, Γ_φu_to_φdeν, Γ_φu_to_φdμν, Γ_φd_to_χq, 
	pi2σ_φuφd_to_quqd, pi2σ_φuφd_to_lvl, pi2σ_φuφd_to_Wpγ, pi2σ_φuφd_to_WpZ, 
	pi2σ_φuφd_to_Wpg, pi2σ_φuφd_to_Wph, pi2σ_φuφu_to_uu, pi2σ_φuφu_to_dd, 
	pi2σ_φuφu_to_ll, pi2σ_φuφu_to_vv, pi2σ_φuφu_to_WW, pi2σ_φuφu_to_γγ, 
	pi2σ_φuφu_to_ZZ, pi2σ_φuφu_to_Zγ, pi2σ_φuφu_to_Zh, pi2σ_φuφu_to_γh, 
	pi2σ_φuφu_to_gg, pi2σ_φuφu_to_gγ, pi2σ_φuφu_to_gZ, pi2σ_φuφu_to_gh,
	pi2σ_φdφd_to_uu, pi2σ_φdφd_to_dd, pi2σ_φdφd_to_ll, pi2σ_φdφd_to_vv, 
	pi2σ_φdφd_to_WW, pi2σ_φdφd_to_γγ, pi2σ_φdφd_to_ZZ, pi2σ_φdφd_to_Zγ, 
	pi2σ_φdφd_to_Zh, pi2σ_φdφd_to_γh, pi2σ_φdφd_to_gg, pi2σ_φdφd_to_gγ, 
	pi2σ_φdφd_to_gZ, pi2σ_φdφd_to_gh


using QuadGK
import PhysicalConstants.CODATA2018: FineStructureConstant

include((@__DIR__) * "../../muc_venv/hubble.jl")
using .Hubble

const αEM = float(FineStructureConstant)
const sw = sqrt(.23121) # sine of weak-mixing angle squared
const cw = sqrt(1 - sw^2)
const mZ = 91.1876 # GeV
const mW = 80.370 # GeV
const mμ = .1056583755 # GeV
const mt = 173.0 # GeV
const α2 = αEM/sw^2
const αS = 0.118
const mh = 125.0
const vev = 246.22

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
Struct that includes relevant parameters for the one-flavor model.
Most observables take this as an input.

λ : Yukawa coupling of the operator φ Q χ

mφ : Mass of the new charged scalar φ

gφ : spin degrees of freedom for φ, should be 1

mχ : Mass of new fermion χ

gχ : spin degrees of freedom for χ, should be 2

h : Hubble evaluated at mφ
"""
struct QLFDM{T <: Real}
	λ::T
	mφ::T 
	gφ::T
	mχ::T
	gχ::T
	δm::T
	h::T
end

function QLFDM(λ::Real, mφ::Real, mχ::Real)
	hub = hubble(mφ)
	dm = real(deltaM(complex(mφ), 2. / 3, -1. / 3, 1. / 6))
	QLFDM(promote(λ, mφ, 1, mχ, 2, dm, hub)...)
end

function QLFDM(λ::Real, mφ::Real, gφ::Real, mχ::Real, gχ::Real)
	hub = hubble(mφ)
	dm = real(deltaM(complex(mφ), 2. / 3, -1. / 3, 1. / 6))
	QLFDM(promote(λ, mφ, gφ, mχ, gχ, dm, hub)...)
end

# Decays 
function Γ_φu_to_χq(model::QLFDM)
	model.λ^2*model.mφ*(1-model.mχ^2/model.mφ^2)^2/(16*π)
end

function Γ_φu_to_χt(model::QLFDM)
	if model.mφ - model.mχ < mt
		return 0.
	else
		sq_factor = sqrt(mt^4 + (model.mφ^2 - model.mχ^2)^2 - 2*mt^2*(model.mφ^2 + model.mχ^2))
		return model.λ^2*(model.mφ^2 - model.mχ^2 - mt^2)*sq_factor / (16*π*model.mφ^3) 
	end
end

function Γ_φu_to_φdeν(model::QLFDM)
    αEM^2*model.δm^5 / (30*π*sw^4*mW^4)
end

function Γ_φu_to_φdμν(model::QLFDM)
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

function Γ_φd_to_χq(model::QLFDM)
	return Γ_φu_to_χq(model)
end

#############################################
# Annihilations
#############################################

function pi2σ_φuφd_to_quqd(x, model::QLFDM)
	w = mW / model.mφ
	3*π*αEM^2*(x-4)^1.5*sqrt(x) / (32*sw^4*(w^2-x)^2)
end

function pi2σ_φuφd_to_lvl(x, model::QLFDM)
	w = mW / model.mφ
	π*αEM^2*(x-4)^1.5*sqrt(x) / (32*sw^4*(w^2-x)^2)
end

function pi2σ_φuφd_to_Wpγ(x, model::QLFDM)
	w = mW / model.mφ
#	logterm1 = (x+4)*log(0.5*(x-sqrt((x-4)*x)-2))
	logterm1 = (x+4)*log((sqrt(x)-sqrt(x-4))/(sqrt(x)+sqrt(x-4)))
#	logterm2 = -4*(x+1)*log(0.5*(x+sqrt((x-4)*x)-2))
	logterm2 = +4*(x+1)*log((sqrt(x)-sqrt(x-4))/(sqrt(x)+sqrt(x-4)))
	pref = π*αEM^2 / (36*sw^2*x^2.5*(x-w^2))
	return pref*(2*sqrt(x-4)*(w^4*(x-3)+3*w^2*(2-3*x)*x+x^2*(x+25))-(w^2-4)*x^1.5*(logterm1 + logterm2))
end

function pi2σ_φuφd_to_WpZ(x, model::QLFDM)
	w = mW / model.mφ
	z = mZ / model.mφ
	pref = π*αEM^2 / (1728*cw^2*sw^4*x^2.5)
	poly = (2*sqrt((x-4)*(w^4 + (x-z^2)^2 - 2*w^2*(x + z^2)))*
           	(8*w^10*(3*(x-4) - 6*sw^2*(x-4) + 4*sw^4*(x-3)) + 
			w^8*(8*x*(48*(sw^2-1)^2 + (-39 + 78*sw^2 - 44*sw^4)*x) + 
			(3*(-4 + x)*(-31 + 14*x) - 6*sw^2*(-4 + x)*(-31 + 14*x) + sw^4*(372 + x*(-277 + 52*x)))*z^2) + 
			3*(-1 + sw^2)^2*(-4 + x)*z^2*(x - z^2)^2*(x^2 + 10*x*z^2 + z^4) + 
			w^6*(32*(9-18*sw^2+10*sw^4)*x^2*(1+2*x) + 4*x*(24 + 3*(50-19*x)*x + sw^4*(24 + 2*(76-31*x)*x) + 6*sw^2*(-8 + x*(19*x-50)))*z^2 + 
			(-528*(-1 + sw^2)^2 + 8*(57 - 114*sw^2 + 58*sw^4)*x - 117*(-1 + sw^2)^2*x^2)*z^4) + 
			2*w^4*(-4*x^3*(168 + 39*x + 44*sw^4*(4 + x) - 6*sw^2*(56 + 13*x)) + 
			x^2*(3*(-412 + x*(171 + 7*x)) - 6*sw^2*(-412 + x*(171 + 7*x)) + sw^4*(-1236 + x*(529 + 26*x)))*z^2 + 
			x*(-720*(-1 + sw^2)^2 + 8*(54 - 108*sw^2 + 53*sw^4)*x - 45*(-1 + sw^2)^2*x^2)*z^4 + 3*(-1 + sw^2)^2*(52 + x*(-37 + 3*x))*z^6) + 
			w^2*(8*x^4*(4*sw^4*(25 + x) + 3*(32 + x) - 6*sw^2*(32 + x)) - 24*x^3*(12 + 9*x - 6*sw^2*(4 + 3*x) + 2*sw^4*(6 + 5*x))*z^2 + 
			x^2*(3*(x-4)*(52+x) - 6*sw^2*(x-4)*(52+x) + sw^4*(-624 + x*(152+3*x)))*z^4 + 6*(sw^2-1)^2*x*(176 + x*(5*x-52))*z^6 + 
			3*(-1 + sw^2)^2*(-16 + x^2)*z^8)))/((w^3 - w*x)^2*(w^4 + (x - z^2)^2 + w^2*(-2*x + (-2 + x)*z^2))) + 
			(4*x^1.5*((-3 + 2*sw^2)*(2*(-4 + w^2)*(w^2 - x)*(-3*(2 + x) + 2*sw^2*(4 + x)) - 
			(sw^2*w^4 + x*(-3*(12 + x) + sw^2*(40 + x)) + w^2*(2*sw^2*(4 - 7*x) + 3*(-4 + 5*x)))*z^2 + 6*(-1 + sw^2)*(-4 + w^2 + x)*z^4 + 
			3*(-1 + sw^2)*z^6)*log((x^1.5 - sqrt(x)*(w^2 + z^2) - sqrt((-4 + x)*(w^4 + (x - z^2)^2 - 2*w^2*(x + z^2))))/
			(x^1.5 - sqrt(x)*(w^2 + z^2) + sqrt((-4 + x)*(w^4 + (x - z^2)^2 - 2*w^2*(x + z^2))))) + 
			(3 - 4*sw^2)*(2*(w^2-4)*(w^2 - x)*(4*sw^2*(1 + x) - 3*(2 + x)) + 
			(sw^2*w^4 + w^2*(12 + 16*sw^2*(-1 + x) - 15*x) + x*(3*(12 + x) - sw^2*(32 + 5*x)))*z^2 + 6*(sw^2-1)*(-4 + w^2 + x)*z^4 + 
			3*(-1 + sw^2)*z^6)*log((sqrt(x)*(w^2 - x + z^2) - sqrt((-4 + x)*(w^4 + (x - z^2)^2 - 2*w^2*(x + z^2))))/
			(sqrt(x)*(w^2 - x + z^2) + sqrt((-4 + x)*(w^4 + (x - z^2)^2 - 2*w^2*(x + z^2)))))))/((w^2 - x)*(w^2 - x + z^2))
	return pref*poly
end

function pi2σ_φuφd_to_Wpg(x, model::QLFDM)
	w = mW / model.mφ
#	logarg = (x - sqrt((x-4)*x)-2) / (x + sqrt((x-4)*x)-2)
	logfactor = (log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
	pref = π*αEM*αS / (3*sw^2*x^1.5*(x-w^2))
	return pref*(2*sqrt(x-4)*(w^4-3*w^2*x+x*(x+4))-(w^2-4)*(x-2)*sqrt(x)*logfactor)
end

function pi2σ_φuφd_to_Wph(x, model::QLFDM)
	w = mW / model.mφ
	xh = mh / model.mφ
	xv = vev / model.mφ
	pref = π^2*αEM^3*xv^2*(x-4)^1.5 / (96*sw^6*w^2*x^2.5*(w^2-x)^2)
	return pref*sqrt(x-(w-xh)^2)*sqrt(x-(w+xh)^2)*(w^4 - 2*xh^2*(w^2+x)+10*w^2*x+x^2+xh^4)
end

# Ignoring processes like φu d -> W χ, since these are all suppressed by λ, 
#   but maybe we should include them?
#   similarly, processes like φu g -> u χ

function pi2σ_φuφu_to_uu(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π * ((x-4)/x)^1.5 / (31104*cw^4*sw^4*(x - z^2)^2)
	term1 = 9*αEM^2*((80*sw^4 - 144*sw^2 + 81)*x^2 - 32*sw^2*(4*sw^4-13*sw^2+9)*x*z^2 + 512*sw^4*(sw^2-1)^2*z^4)
 	term2 = 5184*αS^2*(sw^2-1)^2*sw^4*(x-z^2)^2
 	return pref*(term1 + term2)
end

function pi2σ_φuφu_to_dd(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π * ((x-4)/x)^1.5 / (3456*cw^4*sw^4*(x-z^2)^2)
	term1 = αEM^2*((104*sw^4-180*sw^2+81)*x^2 - 16*sw^2*(8*sw^4 - 17*sw^2 + 9)*x*z^2 + 128*sw^4*(sw^2-1)^2*z^4)
	term2 = 576*αS^2*(sw^2-1)^2*sw^4*(x-z^2)^2
	return pref*(term1 + term2)
end

function pi2σ_φuφu_to_ll(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π*αEM^2*((x-4)/x)^1.5 / (1152*cw^4*sw^4*(x-z^2)^2)
	return pref*(48*sw^2*(sw^2-1)*x*z^2+(8*sw^4-12*sw^2+9)*x^2+128*sw^4*(sw^2-1)^2*z^4)
end

function pi2σ_φuφu_to_vv(x, model::QLFDM)
	z = mZ / model.mφ
	return π*αEM^2*(3-4*sw^2)^2*(x-4)^1.5 * sqrt(x) / (1152*cw^4*sw^4*(x-z^2)^2)
end

function pi2σ_φuφu_to_WW(x, model::QLFDM)
	z = mZ / model.mφ
	w = mW / model.mφ
	pref = π*αEM^2 / (5184*sw^4*x^3*(x-z^2)^2)
	lognum = x - 2*w^2 - sqrt((x-4)*(x-4*w^2))
	logden = x - 2*w^2 + sqrt((x-4)*(x-4*w^2))
	term1den = w^4*(w^4-4*w^2+x)
	term1 = sqrt((x-4)*(x-4*w^2)) * (108*w^4*x^2*(16*w^6 + x^2*(28 + 3*x) - 4*w^2*x*(8 + 5*x) + 4*w^4*(-16 + x*(3 + x))) + 
				12*w^2*x*(4*sw^2*(-4*w^2 + w^4 + x)*(12*w^4*(-8 + x) + (-4 + x)*x^2 + 4*w^2*x*(-26 + 5*x)) - 
				3*x*(12*w^8 + (-4 + x)*x^2 + 2*w^2*x*(40 + 17*x) + w^6*(-80 + 44*x) + w^4*(128 + (-192 + x)*x)))*z^2 + 
				(-16*sw^4*(4*w^2 - x)*(-4 + x)*(-4*w^2 + w^4 + x)*(12*w^4 + 20*w^2*x + x^2) + 
					24*sw^2*x*(-4*w^2 + w^4 + x)*(24*w^6 - 18*w^2*(-4 + x)*x - (-4 + x)*x^2 + 4*w^4*(-16 + 7*x)) + 
					9*x^2*(60*w^8 + 16*w^2*(-4 + x)*x + (-4 + x)*x^2 + 4*w^6*(-92 + 5*x) + w^4*(512 + (-36 + x)*x)))*z^4)/term1den
	term2 = 36*x*(x-z^2)*(z^2*(4*sw^2*(2*w^2+x-8)*(2*w^2+x)+3*x*(2*w^2-x-8))-6*x*(2*w^4+w^2*(3*x-8)-8*x))*log(lognum / logden)
	return pref*( term1 + term2 )
end

function pi2σ_φuφu_to_ZZ(x, model::QLFDM)
	z = mZ / model.mφ
#	lognum = (sqrt((x-4)*(x-4*z^2))-x+2*z^2)^2
#	logden = (sqrt((x-4)*(x-4*z^2))+x+2*z^2)^2
#	pref = π*αEM^2*(3-4*sw^2)^4 / (7776*cw^4*sw^4*x)
#	term1 = (sqrt((x-4)*(x-4*z^2))*(4*x+5*z^4-24*z^2+16))/(z+z^4-4*z^2)
#	term2 = -((4*x*(z^2-2)+(z^2-4)^2)*log(lognum / logden))/(x-2*z^2)
#	return pref*(term1 + term2)
	pref = π*αEM^2*(3-4*sw^2)^4 / (15552*cw^4*sw^4*x)
	logfactor = log((x-2*z^2-sqrt(x-4)*sqrt(x-4*z^2))^2/(x-2*z^2+sqrt(x-4)*sqrt(x-4*z^2))^2)
	polyterm = (sqrt((x-4)*(x-4*z^2))*(4*x+5*z^4-24*z^2+16))/(x+z^4-4*z^2)
	return pref*(polyterm - (4*x*(z^2-2)+(z^2-4)^2)/(x-2*z^2)*logfactor)
end

function pi2σ_φuφu_to_Zγ(x, model::QLFDM)
	z = mZ / model.mφ
#	logarg = (x - sqrt((x-4)*x)-2) / (x + sqrt((x-4)*x)-2)
	logfactor = (log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
	pref = π*αEM^2*(3-4*sw^2)^2 / (162*cw^2*sw^2*x^1.5*(x-z^2))
#	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4)-(x-2)*sqrt(x)*log(logarg))
	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4)-(x-2)*sqrt(x)*logfactor)
end

function pi2σ_φuφu_to_Zh(x, model::QLFDM)
	z = mZ / model.mφ
	xh = mh / model.mφ
	xv = vev / model.mφ	
	pref = π^2*αEM^3*(3-4*sw^2)^2*xv^2*(x-4)^1.5 / (1728*cw^6*sw^6*x^2.5*z^2*(x-z^2)^2)
	return pref*sqrt(x-(xh-z)^2)*sqrt(x-(xh+z)^2)*(2*z^2*(5*x-xh^2)+(x-xh^2)^2+z^4)
end

function pi2σ_φuφu_to_γγ(x, model::QLFDM)
#	4*π*αEM^2*((4+x)*sqrt((x-4)*x) + 2(x-2)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(27*x^2)
	4*π*αEM^2*((4+x)*sqrt((x-4)*x) + 2(x-2)*(log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x)))))/(27*x^2)
end

function pi2σ_φuφu_to_gg(x, model::QLFDM)
#	π*αS^2*((5*x+62)*sqrt(x*(x-4))+4*(4*x+1)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(6*x^2)
	π*αS^2*((5*x+62)*sqrt(x*(x-4))+4*(4*x+1)*(log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x)))))/(6*x^2)
end

function pi2σ_φuφu_to_gγ(x, model::QLFDM)
#	8*π*αEM*αS*((4+x)*sqrt((x-4)*x) + 2(x-2)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(9*x^2)
	8*π*αEM*αS*((4+x)*sqrt((x-4)*x) + 2(x-2)*(log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x)))))/(9*x^2)
end

function pi2σ_φuφu_to_gZ(x, model::QLFDM)
	z = mZ / model.mφ
#	logarg = (x - sqrt((x-4)*x)-2) / (x + sqrt((x-4)*x)-2)
	logfactor = (log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
	pref = 	π*αEM*αS*(3-4*sw^2)^2 / (54*cw^2*sw^2*x^1.5*(x-z^2))
#	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4) - (x-2)*sqrt(x)*(z^2-4)*log(logarg))
	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4) - (x-2)*sqrt(x)*(z^2-4)*logfactor)
end

####################################
# Same as above, but for φd
####################################

function pi2σ_φdφd_to_uu(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π * ((x-4)/x)^1.5 / (31104*cw^4*sw^4*(x-z^2)^2)
	term1 = 9*αEM^2*((116*sw^4-180*sw^2+81)*x^2 - 16*sw^2*(14*sw^4-23*sw^2+9)*x*z^2 + 128*sw^4*(sw^2-1)^2*z^4)
	term2 = 5184*αS^2*(sw^2-1)^2*sw^4*(x-z^2)^2
	return pref*(term1 + term2)
end

function pi2σ_φdφd_to_dd(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π*((x-4)/x)^1.5 / 3456.
	term1 = 576*αS^2
	term2 = (αEM^2/(sw^4*(sw^2-1)^2*(x-z^2)^2))*((68*sw^4-144*sw^2+81)*x^2 - 8*sw^2*(10*sw^4-19*sw^2+9)*x*z^2+32*sw^4*(sw^2-1)^2*z^4)
	return pref*(term1 + term2)
end

function pi2σ_φdφd_to_ll(x, model::QLFDM)
	z = mZ / model.mφ
	pref = π*αEM^2*((x-4)/x)^1.5 / (1152*cw^4*sw^4*(x-z^2)^2)
	return pref*((20*sw^4 - 24*sw^2 + 9)*x^2 - 24*sw^2*(2*sw^4-3*sw^2+1)*x*z^2 + 32*sw^4*(sw^2-1)^2*z^4)
end

function pi2σ_φdφd_to_vv(x, model::QLFDM)
	z = mZ / model.mφ
	return π*αEM^2*(3-2*sw^2)^2*(x-4)^1.5 * sqrt(x) / (1152*cw^4*sw^4*(x-z^2)^2)
end

function pi2σ_φdφd_to_WW(x, model::QLFDM)
	z = mZ / model.mφ
	w = mW / model.mφ
	pref = π*αEM^2 / (5184*sw^4*x^3*(x-z^2)^2)
	lognum = x - 2*w^2 + sqrt((x-4)*(x-4*w^2))
	logden = x - 2*w^2 - sqrt((x-4)*(x-4*w^2))
	term1den = w^4*(w^4-4*w^2+x)
	term1 = sqrt((x-4)*(x-4*w^2)) * (108*w^4*x^2*(16*w^6 + x^2*(28 + 3*x) - 4*w^2*x*(8 + 5*x) + 4*w^4*(-16 + x*(3 + x))) + 
		12*w^2*x*(2*sw^2*(-4*w^2 + w^4 + x)*(12*w^4*(-8 + x) + (-4 + x)*x^2 + 4*w^2*x*(-26 + 5*x)) - 
			3*x*(12*w^8 + (-4 + x)*x^2 + 2*w^2*x*(40 + 17*x) + w^6*(-80 + 44*x) + w^4*(128 + (-192 + x)*x)))*z^2 + 
		(-4*sw^4*(4*w^2 - x)*(-4 + x)*(-4*w^2 + w^4 + x)*(12*w^4 + 20*w^2*x + x^2) + 
			12*sw^2*x*(-4*w^2 + w^4 + x)*(24*w^6 - 18*w^2*(-4 + x)*x - (-4 + x)*x^2 + 4*w^4*(-16 + 7*x)) + 
			9*x^2*(60*w^8 + 16*w^2*(-4 + x)*x + (-4 + x)*x^2 + 4*w^6*(-92 + 5*x) + w^4*(512 + (-36 + x)*x)))*z^4)/term1den
	term2 = 36*x*(x-z^2)*(6*x*(2*w^4 - 8*x + w^2*(-8 + 3*x)) + (3*x*(8 - 2*w^2 + x) - 2*sw^2*(-8 + 2*w^2 + x)*(2*w^2 + x))*z^2)*log(lognum / logden)
	return pref*( term1 + term2 )
end

function pi2σ_φdφd_to_ZZ(x, model::QLFDM)
	z = mZ / model.mφ
#	lognum = (sqrt((x-4)*(x-4*z^2))-x+2*z^2)^2
#	logden = (sqrt((x-4)*(x-4*z^2))+x+2*z^2)^2
#	pref = π*αEM^2*(3-2*sw^2)^4 / (7776*cw^4*sw^4*x)
#	term1 = (sqrt((x-4)*(x-4*z^2))*(4*x+5*z^4-24*z^2+16))/(z+z^4-4*z^2)
#	term2 = -((4*x*(z^2-2)+(z^2-4)^2)*log(lognum / logden))/(x-2*z^2)
#	return pref*(term1 + term2)
	pref = π*αEM^2*(3-2*sw^2)^4 / (15552*cw^4*sw^4*x)
	logfactor = log((x-2*z^2-sqrt(x-4)*sqrt(x-4*z^2))^2/(x-2*z^2+sqrt(x-4)*sqrt(x-4*z^2))^2)
	polyterm = (sqrt((x-4)*(x-4*z^2))*(4*x+5*z^4-24*z^2+16))/(x+z^4-4*z^2)
	return pref*(polyterm - (4*x*(z^2-2)+(z^2-4)^2)/(x-2*z^2)*logfactor)
end

function pi2σ_φdφd_to_Zγ(x, model::QLFDM)
	z = mZ / model.mφ
#	logarg = (x - sqrt((x-4)*x)-2) / (x + sqrt((x-4)*x)-2)
	logfactor = (log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
	pref = π*αEM^2*(3-2*sw^2)^2 / (648*cw^2*sw^2*x^1.5*(x-z^2))
#	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4)-(x-2)*sqrt(x)*log(logarg))
	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4)-(x-2)*sqrt(x)*logfactor)
end

function pi2σ_φdφd_to_Zh(x, model::QLFDM)
	z = mZ / model.mφ
	xh = mh / model.mφ
	xv = vev / model.mφ
	pref = π^2*αEM^3*(3-2*sw^2)^2*xv^2*(x-4)^1.5 / (1728*cw^6*sw^6*x^2.5*z^2*(x-z^2)^2)
	return pref*sqrt(x-(xh-z)^2)*sqrt(x-(xh+z)^2)*(2*z^2*(5*x-xh^2)+(x-xh^2)^2+z^4)
end

function pi2σ_φdφd_to_γγ(x, model::QLFDM)
#	π*αEM^2*((4+x)*sqrt((x-4)*x) + 2(x-2)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(108*x^2)
	π*αEM^2*((4+x)*sqrt((x-4)*x) + 2(x-2)*(log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x)))))/(108*x^2)
end

function pi2σ_φdφd_to_gg(x, model::QLFDM)
#	π*αS^2*((5*x+62)*sqrt(x*(x-4))+4*(4*x+1)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(6*x^2)
	π*αS^2*((5*x+62)*sqrt(x*(x-4))+4*(4*x+1)*(log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x)))))/(6*x^2)	
end

function pi2σ_φdφd_to_gγ(x, model::QLFDM)
#	2*π*αEM*αS*((4+x)*sqrt((x-4)*x) + 2(x-2)*log((x-sqrt(x*(x-4))-2)/(x+sqrt(x*(x-4))-2)))/(9*x^2)
	2*π*αEM*αS*(sqrt(x*(x-4))*(x+4)+2*(x-2)*log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - 2*(x-2)*log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))/(9*x^2)
end

function pi2σ_φdφd_to_gZ(x, model::QLFDM)
	z = mZ / model.mφ
#	logarg = (x - sqrt((x-4)*x)-2) / (x + sqrt((x-4)*x)-2)
	logfactor = (log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
	pref = 	π*αEM*αS*(3-2*sw^2)^2 / (54*cw^2*sw^2*x^1.5*(x-z^2))
#	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4) - (x-2)*sqrt(x)*(z^2-4)*log(logarg))
	return pref*(2*sqrt(x-4)*(-3*x*z^2+x*(x+4)+z^4) - (x-2)*sqrt(x)*(z^2-4)*logfactor)
end

end