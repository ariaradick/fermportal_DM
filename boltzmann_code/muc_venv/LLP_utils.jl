module LLP_Utils

export bin_data, RunsInfo, RunSummary, N_events, mass, initial_state, beta,
        betagamma, eta, β, βγ, η, xsec, xsec_nocoef, Lτ, Lzτ, dir_runx,
        get_runs_info, get_run_summaries, eid_mchi, eid_tauphi, nopdf_run_info

using NaturallyUnitful
using CSV
using DataFrames
using Interpolations

RUN_RESULTS_FILENAME = "run_info.csv"
SUMMARY_FILENAME = "diphi_summary.csv"

to_cm_ns = ustrip(unnatural(u"cm*ns^-1", 1))
GeV_to_ns = ustrip(unnatural(u"ns", 1*u"GeV^-1"))
ccms = ustrip(unnatural(u"cm*s^-1", 1))

#~~~ Dynamic Coefficients ~~~
m_mu = .10566
cqgrid = CSV.read((@__DIR__) * "/c_vs_Q_grid_resc1em4.csv", DataFrame)
tt = 2*log.(cqgrid.Q/m_mu)
itp_cL = interpolate((tt,), cqgrid.cL, Gridded(Linear()))
itp_cR = interpolate((tt,), cqgrid.cR, Gridded(Linear()))

f_cL(Q) = itp_cL(2*log(Q/m_mu))
f_cR(Q) = itp_cR(2*log(Q/m_mu))

coefLR(Q) = f_cL(Q)*f_cR(Q)
coefL(Q) = f_cL(Q)
coefR(Q) = f_cR(Q)
coef1(Q) = 1.0

coef_fn_LR = [coefLR, coefR, coefL, coef1]
coef_fn_RL = [coefLR, coefL, coefR, coef1]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sqrt(s)/2 coefficients
coefs = Dict("1" => 1.0, "L" => 0.2862, "R" => 0.2673, "LR" => 0.2673*0.2862)
coefs_LR = [coefs["LR"], coefs["R"], coefs["L"], coefs["1"]]
coefs_RL = [coefs["LR"], coefs["L"], coefs["R"], coefs["1"]]

function bin_data(x, binedges)
    N_bins = length(binedges)-1
    res = zeros(eltype(x), N_bins)
    for i in 1:N_bins
        res[i] = count(y -> binedges[i] <= y < binedges[i+1], x)
    end
    res[end] = count(y -> binedges[end-1] <= y <= binedges[end], x)
    return res
end

struct RunsInfo
    mphis::Vector{Float64}
    init_states::Vector{String}
    dirs::Matrix{String}
    parent_dir::String
    xsecs::Matrix{Float64}
    coefs::Matrix{Float64}
end

struct RunSummary
    β::Vector{Float64}
    βγ::Vector{Float64}
    η::Vector{Float64}
    xsec::Float64
    initial_state::String
    mphi::Float64
    coef::Float64
end

function N_events(run_summary::RunSummary)
    length(run_summary.β)
end

function mass(run_summary::RunSummary)
    run_summary.mphi
end

function initial_state(run_summary::RunSummary)
    run_summary.initial_state
end

function beta(run_summary::RunSummary)
    run_summary.β
end

β(x::RunSummary) = beta(x)

function betagamma(run_summary::RunSummary)
    run_summary.βγ
end

βγ(x::RunSummary) = betagamma(x)

function eta(run_summary::RunSummary)
    run_summary.η
end

η(x::RunSummary) = eta(x)

function xsec(run_summary::RunSummary)
    run_summary.xsec
end

function xsec_nocoef(run_summary::RunSummary)
    run_summary.xsec / run_summary.coef
end

function Lτ(run_summary::RunSummary)
    βγ(run_summary) ./ cosh.(η(run_summary))
end

function Lzτ(run_summary::RunSummary)
    βγ(run_summary) .* tanh.(η(run_summary))
end

function dir_runx(x)
    if x < 10
        return "Events/run_0$x/"
    else
        return "Events/run_$x/"
    end
end

function load_summary(rundir, xsec, init_state, mphi, coef; sfname=SUMMARY_FILENAME)
    df = CSV.read(rundir*sfname, DataFrame)
    return RunSummary(Vector(df.beta), Vector(df.betagamma), Vector(df.eta),
                      xsec, init_state, mphi, coef)
end

function load_run_info(mgdir)
    CSV.read(mgdir*RUN_RESULTS_FILENAME, DataFrame)
end

function get_runs_info(parent_dir, LR_dir, RL_dir, vbf_dir; dyn=true, 
                       vv_rts=false, sfname=SUMMARY_FILENAME)

    LR_info = load_run_info(parent_dir*LR_dir)
    RL_info = load_run_info(parent_dir*RL_dir)
    vbf_info = load_run_info(parent_dir*vbf_dir)

    lpp_pairs = [(0,0), (0,1), (1,0), (1,1)]
    LR_idxs = Vector{Vector{Int}}()
    RL_idxs = Vector{Vector{Int}}()
    for lpp in lpp_pairs
        push!(LR_idxs, findall(x -> x == "$lpp", LR_info[:,1]))
        push!(RL_idxs, findall(x -> x == "$lpp", LR_info[:,1]))
    end

    if vv_rts
        vbf_idxs = findall(x -> x == 10e3, vbf_info[:,1]) # sqrt(s) for the vbf runs
    else
        vbf_idxs = 1:length(vbf_info.xsec)
    end

    LR_infos = [LR_info[idxs,:] for idxs in LR_idxs]
    LR_mphis = [minfo.mphi for minfo in LR_infos]
    if dyn
        LR_coefs = [coef_fn_LR[i].(LR_mphis[i]) for i in eachindex(LR_infos)]
    else
        LR_coefs = [coefs_LR[i].*ones(Float64, size(LR_mphis[i])) for i in eachindex(LR_infos)]
    end
    LR_xsecs = [LR_infos[i].xsec .* LR_coefs[i] for i in eachindex(LR_infos)]

    RL_infos = [RL_info[idxs,:] for idxs in RL_idxs]
    RL_mphis = [minfo.mphi for minfo in RL_infos]
    if dyn
        RL_coefs = [coef_fn_RL[i].(RL_mphis[i]) for i in eachindex(RL_infos)]
    else
        RL_coefs = [coefs_RL[i].*ones(Float64, size(RL_mphis[i])) for i in eachindex(RL_infos)]
    end
    RL_xsecs = [RL_infos[i].xsec .* RL_coefs[i] for i in eachindex(RL_infos)]

    vbf_info = vbf_info[vbf_idxs,:]

    LR_init = ["LR$(lpp...)" for lpp in lpp_pairs]
    RL_init = ["RL$(lpp...)" for lpp in lpp_pairs]
    init_states = [LR_init..., RL_init..., "VV"]

    all_ids = [LR_idxs..., RL_idxs..., vbf_idxs]
    topdirs = [LR_dir, LR_dir, LR_dir, LR_dir, RL_dir, RL_dir, RL_dir, RL_dir,
               vbf_dir]
    
    all_xsecs = [LR_xsecs..., RL_xsecs..., vbf_info.xsec]
    all_coefs = [LR_coefs..., RL_coefs..., ones(size(vbf_info.mphi))]

    run_dirs = Matrix{String}(undef, length(all_ids), length(vbf_info.mphi))
    xsecs = zeros(Float64, (length(all_ids), length(vbf_info.mphi)))
    coefs = zeros(Float64, (length(all_ids), length(vbf_info.mphi)))
    for i in eachindex(all_ids)
        run_dirs[i,:] = topdirs[i].*dir_runx.(all_ids[i])
        xsecs[i,:] = all_xsecs[i]
        coefs[i,:] = all_coefs[i]
    end

    for m in LR_mphis
        if m ≉ vbf_info.mphi
            error("mphis don't match!")
        end
    end

    for m in RL_mphis
        if m ≉ vbf_info.mphi
            error("mphis don't match!")
        end
    end

    return RunsInfo(vbf_info.mphi, init_states, run_dirs, parent_dir, xsecs,
                    coefs)
end

function get_run_summaries(run_infos::RunsInfo, mass_idx; sfname=SUMMARY_FILENAME)
    load_summary.(run_infos.parent_dir.*run_infos.dirs[:,mass_idx], 
                  run_infos.xsecs[:,mass_idx], run_infos.init_states,
                  run_infos.mphis[mass_idx]', run_infos.coefs[:,mass_idx];
                  sfname=sfname)
end

function get_run_summaries(run_infos::RunsInfo; sfname=SUMMARY_FILENAME)
    load_summary.(run_infos.parent_dir.*run_infos.dirs, 
                  run_infos.xsecs, run_infos.init_states, run_infos.mphis',
                  run_infos.coefs; sfname=sfname)
end

function nopdf_run_info(run_infos::RunsInfo; idxs=[1,5,9])
    new_cs = ones(Float64, (3,length(run_infos.mphis)))
    new_cs[1:2,:] = fill(0.25, (2,length(run_infos.mphis)))

    RunsInfo(run_infos.mphis, run_infos.init_states[idxs], 
    run_infos.dirs[idxs,:], run_infos.parent_dir,
    new_cs.*run_infos.xsecs[idxs,:]./run_infos.coefs[idxs,:], new_cs)
end

function Nexp_decays(summary::RunSummary, lumi, tauphi, R_cm, Z; Z_is_eta=false)
    pref = summary.xsec * lumi / N_events(summary)
    L = Lτ(summary).*tauphi.*to_cm_ns
    if !(Z_is_eta)
        z = Lzτ(summary).*tauphi.*to_cm_ns
    else
        z = eta(summary)
    end
    return pref*count((R_cm[1] .< L .< R_cm[2]) .&& (z .< Z))
end

detector_Rs = [[3, 10.4], [12.7, 55.4], [81.9, 148.6], [150, 170.2], 
                    [174, 333], [348.3, 429.0], [446.1, 645.0], [645.0, Inf]]
detector_Zs = [65, (69.2+48.2)/2, 124.9, 221.0, 221.0, 412.9, 417.9, Inf]

det_labs = ["in Vertex Detector", "in Inner Tracker", "in Outer Tracker", "in ECAL",
            "in HCAL", "in Solenoid", "in Muon Detector", "Outside"]

function _events_in_detector(run_info::RunsInfo, lumi, τφs; 
                            det_comps=(1:length(detector_Rs)))

    summaries = get_run_summaries(run_info)
    ndet = length(det_comps)
    nx = size(τφs)[2]
    nmphi = length(run_info.mphis)
    ninit = length(run_info.init_states)

    res = zeros(Float64, (ndet, nx, ninit, nmphi))
    for (i,idx) in enumerate(det_comps)
        if idx==8
            use_eta=true
            z = 0.61
        else
            use_eta = false
            z = detector_Zs[idx]
        end
        for j in 1:nx
            res[i,j,:,:] = Nexp_decays.(summaries, lumi, τφs[:,j]', 
                                        (detector_Rs[idx],), z; Z_is_eta=use_eta)
        end
    end
    return res
end

function eid_mchi(run_info::RunsInfo, lumi, mchis, tauphis_ns; 
                  det_comps=(1:length(detector_Rs)))
    τφs = tauphis_ns # ns
    eid = _events_in_detector(run_info, lumi, τφs; det_comps=det_comps)
    return (N_exp=eid, rts=10e3, lumi=lumi, mphis=run_info.mphis, mchis=mchis,
            labels=det_labs[det_comps])
end

function eid_tauphi(run_info::RunsInfo, lumi, tauphis_ns; 
                    det_comps=(1:length(detector_Rs)))
    τφs = transpose(reduce(hcat, fill(tauphis_ns, length(run_info.mphis))))
    eid = _events_in_detector(run_info, lumi, τφs; det_comps=det_comps)
    return (N_exp=eid, rts=10e3, lumi=lumi, mphis=run_info.mphis, tauphis=τφs,
            labels=det_labs[det_comps])
end

end
