module Prompt

export BkgSummary, RunsInfo, RunSummary, N_events, mass, initial_state,
        mll, ptll, mt2, xsec, load_summary, load_bkg_summary, load_run_info,
        get_runs_info, get_run_summaries, get_all_run_summaries, fcut,
        N_pass, significances, bin_data

using CSV
using DataFrames

# filenames from summarize_runs.jl
RUN_RESULTS_FILENAME = "run_info.csv"
SUMMARY_FILENAME = "dilepton_summary.csv"

function dir_runx(x)
    if x < 10
        return "Events/run_0$x/"
    else
        return "Events/run_$x/"
    end
end

abstract type Summary end

struct BkgSummary <: Summary
    mll::Vector{Float64}
    ptll::Vector{Float64}
    mt2::Vector{Float64}
    xsec::Float64
end

struct RunsInfo
    mphis::Vector{Float64}
    init_states::Vector{String}
    dirs::Matrix{String}
    parent_dir::String
    xsecs::Matrix{Float64}
    bkg::BkgSummary
end

struct RunSummary <: Summary
    mll::Vector{Float64}
    ptll::Vector{Float64}
    mt2::Vector{Float64}
    xsec::Float64
    initial_state::String
    mphi::Float64
end

#~~~ convenience functions: ~~~
function N_events(run_summary::Summary)
    length(run_summary.mll)
end

function mass(run_summary::RunSummary)
    run_summary.mphi
end

function initial_state(run_summary::RunSummary)
    run_summary.initial_state
end

function mll(run_summary::Summary)
    run_summary.mll
end

function ptll(run_summary::Summary)
    run_summary.ptll
end

function mt2(run_summary::Summary)
    run_summary.mt2
end

function xsec(run_summary::Summary)
    run_summary.xsec
end
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function load_summary(rundir, xsec, init_state, mphi)
    df = CSV.read(rundir*SUMMARY_FILENAME, DataFrame)
    return RunSummary(Vector(df.m_ll), Vector(df.pt_ll), Vector(df.MT2),
                      xsec, init_state, mphi)
end

function load_bkg_summary(rundir, xsec)
    df = CSV.read(rundir*SUMMARY_FILENAME, DataFrame)
    return BkgSummary(Vector(df.m_ll), Vector(df.pt_ll), Vector(df.MT2), xsec)
end

function load_run_info(mgdir)
    CSV.read(mgdir*RUN_RESULTS_FILENAME, DataFrame)
end

# factor of 2 because we're not using RL:
coefs_LR = 2 .* [0.00606, 0.055, 0.055, 0.5]

function get_runs_info(parent_dir, bkg_dir, LR_dir, RL_dir, vbf_dir)

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

    LR_infos = [LR_info[idxs,:] for idxs in LR_idxs]
    LR_xsecs = [LR_infos[i].xsec*coefs_LR[i] for i in eachindex(LR_infos)]
    LR_mphis = [minfo.mphi for minfo in LR_infos]

    RL_infos = [RL_info[idxs,:] for idxs in RL_idxs]
    RL_xsecs = [RL_infos[i].xsec*coefs_RL[i] for i in eachindex(RL_infos)]
    RL_mphis = [minfo.mphi for minfo in RL_infos]

    LR_init = ["LR$(lpp...)" for lpp in lpp_pairs]
    RL_init = ["RL$(lpp...)" for lpp in lpp_pairs]
    init_states = [LR_init..., RL_init..., "VV"]

    all_ids = [LR_idxs..., RL_idxs..., vbf_idxs]
    topdirs = [LR_dir, LR_dir, LR_dir, LR_dir, RL_dir, RL_dir, RL_dir, RL_dir,
               vbf_dir]
    all_xsecs = [LR_xsecs..., RL_xsecs..., vbf_info.xsec]

    run_dirs = Matrix{String}(undef, length(all_ids), length(vbf_info.mphi))
    xsecs = zeros(Float64, (length(all_ids), length(vbf_info.mphi)))
    for i in eachindex(all_ids)
        run_dirs[i,:] = topdirs[i].*dir_runx.(all_ids[i])
        xsecs[i,:] = all_xsecs[i]
    end

    bkg_info = load_run_info(parent_dir*bkg_dir)
    bkg_summary = load_bkg_summary(parent_dir*bkg_dir*dir_runx(1),
                                   bkg_info.xsec[1])

    return RunsInfo(vbf_info.mphi, init_states, run_dirs, parent_dir, xsecs, 
                    bkg_summary)
end

function get_run_summaries(run_infos::RunsInfo, mass_idx)
    load_summary.(run_infos.parent_dir.*run_infos.dirs[:,mass_idx], 
                  run_infos.xsecs[:,mass_idx], run_infos.init_states,
                  run_infos.mphis[mass_idx]')
end

function get_all_run_summaries(run_infos::RunsInfo)
    load_summary.(run_infos.parent_dir.*run_infos.dirs, 
                  run_infos.xsecs, run_infos.init_states, run_infos.mphis')
end

function bin_data(x, binedges)
    N_bins = length(binedges)-1
    res = zeros(eltype(x), N_bins)
    for i in 1:N_bins
        res[i] = count(y -> binedges[i] <= y < binedges[i+1], x)
    end
    res[end] = count(y -> binedges[end-1] <= y <= binedges[end], x)
    return res
end

gtr_min(summ_vec, min) = summ_vec .> min

function fcut(summary::Summary, mll_min, ptll_min, mt2_min)
    N_mll = length(mll_min)
    N_ptll = length(ptll_min)
    N_mt2 = length(mt2_min)
    Nevent = N_events(summary)

    res = zeros(Float64, (N_mll, N_ptll, N_mt2))

    mll_test = gtr_min.( (summary.mll,), mll_min )
    ptll_test = gtr_min.( (summary.ptll,), ptll_min )
    mt2_test = gtr_min.( (summary.mt2,), mt2_min )

    @Threads.threads for i in 1:N_mll
        @Threads.threads for j in 1:N_ptll
            @Threads.threads for k in 1:N_mt2
                res[i,j,k] = count(mll_test[i] .&& ptll_test[j] .&& mt2_test[k])
            end
        end
    end

    return res ./ Nevent
end

fcut(summary::Summary, mins) = fcut(summary, mins...)[1]

function N_pass(summary::Summary, lumi, mins)
    summary.xsec * lumi * fcut(summary, mins)
end

function significances(run_info::RunsInfo, lumi, mll_min, ptll_min, mt2_min)
    mins = (mll_min, ptll_min, mt2_min)

    N_bkg_true = N_pass(run_info.bkg, lumi, mins)
    N_bkg = N_bkg_true .+ ((N_bkg_true .<= 2.0) .* (2.0 .- N_bkg_true))

    summaries = get_all_run_summaries(run_info)

    N_each = N_pass.(summaries, lumi, (mins,))
    N_mass = sum(N_each; dims=1)
    @. N_mass *= (N_mass >= 20)

    return N_mass ./ sqrt(N_bkg)
end

end