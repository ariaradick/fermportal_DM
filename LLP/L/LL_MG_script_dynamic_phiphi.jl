import DelimitedFiles: readdlm

MADGRAPH_EXE = "/home/aradick/Downloads/Madgraph/3_5_6/MG5_aMC_v3_5_6/bin/mg5_aMC"
MADGRAPH_OUTPUT = "/home/aradick/Downloads/Madgraph/3_5_6/"

Mphi = readdlm((@__DIR__)*"/"*"mphi_list.csv")
Mchi = 1e-6 # GeV
Lambda = 1e-7

Ebeam = 5e3 # for sqrt(s) = 10 TeV

# for muon initial states
Hels = ["LR", "RL"]
LHAIDs = [40000204, 40000304]
LPPs = [(0,0), (0,1), (1,0), (1,1)]

DIR_φφ = "mumu_to_phiphi"
DIR_VVφφ = "VV_to_phiphi"

FOLDER = "LL_dynamic/"

mumu_script_file = DIR_φφ * "_script.txt"
vv_script_file = DIR_VVφφ * "_script.txt"

function mumu_cards(dir, lam, mphi, mchi, lpp1, lpp2, lhaid, fact_scale;
               sde_strategy=2, nevents=10000, ebeam=Ebeam)
    """launch $dir
    set mchi $mchi
    set msca $mphi
    set lambda1 $lam
    set lambda2 $lam
    set lambda3 $lam
    set nevents $nevents
    set lpp1 $lpp1
    set lpp2 $lpp2
    set pdlabel lhapdf
    set lhaid $lhaid
    set ebeam1 $ebeam
    set ebeam2 $ebeam
    set fixed_fac_scale True
    set dsqrt_q2fact1 $fact_scale
    set dsqrt_q2fact2 $fact_scale
    set sde_strategy $sde_strategy
    0
    """
end

function mumu_script(mphis, lam, mchi, hels, lhaids, lpp_pairs; 
                     sde_strategy=2, nevents=10000, dir=DIR_φφ)
    res = ""
    for hi in eachindex(hels)
        hp = hels[hi][1]
        hm = hels[hi][2]

        hdir = "$(dir)_$(hels[hi])"

        res *= """import model L_UFO
            generate mu+{$hp} mu-{$hm} > scap scam QED<=3
            output $hdir
            """
        
        for lpps in lpp_pairs
            for mphi in mphis
                res *= mumu_cards(hdir, lam, mphi, mchi, lpps[1], lpps[2], 
                       lhaids[hi], mphi; sde_strategy=sde_strategy,
                       nevents=nevents)
            end
        end

        res *= """launch $hdir -i
        print_results --path=./$hdir/xsec.txt --format=short
        exit
        """
    end

    res *= "exit"

    return res
end

function vv_cards(dir, lam, mphi, mchi; sde_strategy=2, nevents=10000, ebeam=Ebeam)
    """launch $dir
    set mchi $mchi
    set msca $mphi
    set lambda1 $lam
    set lambda2 $lam
    set lambda3 $lam
    set nevents $nevents
    set lpp1 -4
    set lpp2 4
    set fixed_fac_scale false
    set dynamical_scale_choice 4
    set scalefact 0.5
    set ievo_eva 0
    set ebeam1 $ebeam
    set ebeam2 $ebeam
    set no_parton_cut
    set dsqrt_shat 1000
    set sde_strategy $sde_strategy
    0
    """
end

function vv_script(mphis, lam, mchi; sde_strategy=2, nevents=10000, dir=DIR_VVφφ)
    res = """import model LFDM_LL_UFO
    set group_subprocesses false
    define vxp = w+ z a
    define vxm = w- z a
    generate vxp vxm > scap scam QED<=3 LFDM=0
    add process vxp vxm > scap sca0 QED<=3 LFDM=0
    add process vxp vxm > scam sca0~ QED<=3 LFDM=0 
    output $dir
    """

    for mphi in mphis
        res *= vv_cards(dir, lam, mphi, mchi; sde_strategy=sde_strategy,
                       nevents=nevents)
    end

    res *= """launch $dir -i
    print_results --path=./$dir/xsec.txt --format=short
    exit
    exit"""

    return res
end

function main()
    output = MADGRAPH_OUTPUT * FOLDER
    if isdir(output)
        rm(output; recursive=true)
    end
    mkdir(output)

    write(output*mumu_script_file, 
          mumu_script(Mphi, Lambda, Mchi, Hels, LHAIDs, LPPs))
    write(output*vv_script_file, vv_script(Mphi, Lambda, Mchi))

    run(Cmd(`$MADGRAPH_EXE $mumu_script_file`, dir=output))
    run(Cmd(`$MADGRAPH_EXE $vv_script_file`, dir=output))
end

main()
