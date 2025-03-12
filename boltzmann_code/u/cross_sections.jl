module uRFlavorDM

export uRFDM, Γ_φ_to_χq, pi2σ_φφ_to_γγ, pi2σ_φφ_to_γZ, pi2σ_φφ_to_ZZ, 
	pi2σ_φφ_to_WW, pi2σ_φφ_to_gg, pi2σ_φφ_to_gγ, pi2σ_φφ_to_gZ, 
	pi2σ_φφ_to_qq, pi2σ_φφ_to_ll

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
const α2 = αEM/sw^2
const αS = 0.118

function arctanh(x)
    .5*(log(1+x) - log(1-x))
end

function arccoth(x)
    .5*(log(x+1) - log(x-1))
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
struct uRFDM{T <: Real}
	λ::T
	mφ::T
	gφ::T
	mχ::T
	gχ::T
	h::T
end

function uRFDM(λ::Real, mφ::Real, mχ::Real)
	hub = hubble(mφ)
	uRFDM(promote(λ, mφ, 1, mχ, 2, hub)...)
end

function uRFDM(λ::Real, mφ::Real, gφ::Real, mχ::Real, gχ::Real)
	hub = hubble(mφ)
	uRFDM(promote(λ, mφ, gφ, mχ, gχ, hub)...)
end

function Γ_φ_to_χq(model::uRFDM)
	model.λ^2*(model.mφ^2 - model.mχ^2)*(model.mφ^2 + model.mχ^2)/(16*π*model.mφ^3)
end

function pi2σ_φφ_to_γγ(x, model::uRFDM; approx=true)
	4*π*αEM^2 / (27*x^2)*(sqrt(x*(x-4))*(x+4)+2*(x-2)*log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - 2*(x-2)*log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
end

function pi2σ_φφ_to_γZ(x, model::uRFDM; approx=true)
	z = mZ / model.mφ
	pref = 8*π*sw^2*αEM^2/(81*cw^2*x^1.5*(x-z^2))
	t1 = 2*sqrt(x-4)*(x^2 + z^4 + x*(4 - 3*z^2))
	t2 = -(x-2)*sqrt(x)*(z^2-4)*(log((-sqrt(x-4) + sqrt(x))/(sqrt(x-4) + sqrt(x))) - log(-((sqrt(x-4) + sqrt(x))/(sqrt(x-4) - sqrt(x)))))
	return pref*(t1 + t2)
end

function pi2σ_φφ_to_ZZ(x, model::uRFDM; approx=true)
	z = mZ / model.mφ
	kinef = sqrt( (x-4) * (x-4*z^2) )
	pref = 8*π*αEM^2*sw^4/(243*cw^4*x)
	t1 = (kinef*(16+4*x-24*z^2+5*z^4))/(x-4*z^2+z^4) 
	t2 = -(((z^2-4)^2 + 4*x*(z^2-2))*(log((x-2*z^2-kinef)/(x-2*z^2+kinef)) - log((x-2*z^2+kinef)/(x-2*z^2-kinef))))/(x-2*z^2)
	return pref*(t1 + t2)
end

function pi2σ_φφ_to_WW(x, model::uRFDM; approx=true)
	z = mZ / model.mφ
	w = mW / model.mφ

	(αEM^2*π*z^4*(x-4)^(1.5)*sqrt(x-4*w^2)*(x^3 + 16*w^2*x^2 - 68*w^4*x - 48*w^6))/(324*w^4*x^3*(x-z^2)^2)
end

function pi2σ_φφ_to_gg(x, model::uRFDM; approx=true)
    (π*αS^2*(sqrt(x-4)*sqrt(x)*(62+5*x) + 
	  4*(1+4*x)*log((-sqrt(x-4)+sqrt(x))/(sqrt(x-4)+sqrt(x))) - 
	  4*(1+4*x)*log(-((sqrt(x-4)+sqrt(x))/(sqrt(x-4)-sqrt(x))))))/(6*x^2)
end

function pi2σ_φφ_to_gγ(x, model::uRFDM; approx=true)
	8*π*αEM*αS/(9*x^2)*(sqrt(x*(x-4))*(x+4)+2*(x-2)*log((sqrt(x)-sqrt(x-4))/(sqrt(x-4)+sqrt(x))) - 2*(x-2)*log((-sqrt(x-4)-sqrt(x))/(sqrt(x-4)-sqrt(x))))
end

function pi2σ_φφ_to_gZ(x, model::uRFDM; approx=true)
	pi2σ_φφ_to_γZ(x, model::uRFDM, approx=approx)*(3*αS / αEM)
end

function pi2σ_φφ_to_qq(x, model::uRFDM; approx=true)
	z = mZ / model.mφ

	(π*((-4+x)/x)^1.5*(18*(17*x^2+40*(-1+sw^2)*x*z^2 + 32*(-1+sw^2)^2*z^4)*αEM^2 + 648*(-1+sw^2)^2*(x-z^2)^2*αS^2))/(3888*cw^4*(x-z^2)^2)
end

function pi2σ_φφ_to_ll(x, model::uRFDM; approx=true)
	z = mZ / model.mφ
	pref = π*αEM^2/(72*cw^4)
	kinef = ((x-4)/x)^1.5
	return pref*kinef*(8*cw^4*z^4 - 12*cw^2*x*z^2+5*x^2)/((x-z^2)^2)
	#	1 + (π*((x-4)/x)^1.5*(5*x^2+12*(-1+sw^2)*x*z^2 + 8*(-1+sw^2)^2*z^4)*αEM^2)/(72*cw^4*(x-z^2)^2)
end

end
