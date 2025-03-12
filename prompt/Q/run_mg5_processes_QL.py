#!/usr/bin/env python3

#MGDIR="/nfs/user/sdh99/hep-packages/MG5_aMC_v3_5_6"
MGDIR="~/work/MG5_aMC_v3_5_3"

# imports go here
import os, sys
import numpy as np

sys.path.insert(0, 'model_inputs')
from model_inputs_QL import *

resc_lhaid = { "1e-2": 40000202, "1e-4": 40000204, "1e-6": 40000206, "1e-8": 40000208 }

#mphi_values = np.concatenate(( np.arange(100., 350., 25.), np.arange(350., 800., 50.), np.arange(800., 1600., 100.), np.arange(1600., 5000., 200.) )).tolist()

# loading Aria's files for later use
# UPDATE THESE (the phi, chi grids are fine, for now)
mphi_values = np.genfromtxt("../relic_abundance/uR_model/mphis_grid.csv")
mchi_values = np.genfromtxt("../relic_abundance/uR_model/mchis_grid.csv")
#lambda_values = np.genfromtxt("../relic_abundance/QL_model/QL_lambdas_grid.csv")

#mphi_values = np.array([100., 500., 1000., 5000.])

base_script = """set nevents 10000
set ebeam1 5000.0
set ebeam2 5000.0
set use_syst False
set mchi 1e-6
set lambda1 1e-7
set lambda2 1e-7
set lambda3 1e-7
"""

eva_script = """set pdlabel eva
set lpp1 -4
set lpp2 4
set fixed_fac_scale false
set dynamical_scale_choice 4
set scalefact 0.5
set ievo_eva 0
set dsqrt_q2fact1 5000.0
set dsqrt_q2fact2 5000.0
set sde_strategy 1
set nhel 1
set width all 0
"""

def mumu_generate_line(hel1, hel2):
    return """generate mu+{{{0}}} mu-{{{1}}} > scau scaub > j j chi chi~ QED<=3
add process mu+{{{0}}} mu-{{{1}}} > scad scadb > j j chi chi~ QED<=3
""".format(hel1, hel2)

vv_generate_line = """generate vxp vxm > scau scaub > j j chi chi~ QED<=3
add process vxp vxm > scad scadb > j j chi chi~ QED<=3
add process vxp vxm > scau scadb > j j chi chi~ QED<=3
add process vxp vxm > scad scaub > j j chi chi~ QED<=3
"""
 

def create_script_file(proc, hel1, hel2, pdf_resc='1e-4', fact_scale=None):

    filename = ""
    generate_line = ""
    run_name = ""
    extra_commands = ""

    if proc == 'mumu':
        if fact_scale == 'dyn':
            filename = "madgraph_run_files_QL/run_mu{}mu{}_to_phiphibar_chichijj_resc{}_fact_dyn.txt".format(hel1, hel2, pdf_resc)
            generate_line = mumu_generate_line(hel1, hel2)
            run_name = "fermportal_QL_mu{0}mu{1}_to_phiphibar_chichijj_resc{2}_fact_dyn".format(hel1, hel2, pdf_resc)

        else: # the default (factorization scale = beam energy )
            filename = "madgraph_run_files_QL/run_mu{}mu{}_to_phiphibar_chichijj_resc{}.txt".format(hel1, hel2, pdf_resc)
            generate_line = mumu_generate_line(hel1, hel2)
            run_name = "fermportal_QL_mu{0}mu{1}_to_phiphibar_chichijj_resc{}".format(hel1, hel2, pdf_resc)

    elif proc == 'vv':
        filename = 'madgraph_run_files_QL/run_vv_to_phiphibar_chichijj.txt'
        generate_line = vv_generate_line
        run_name = 'fermportal_QL_vv_to_phiphibar_chichijj'

    elif proc == 'aa':
        pass

    elif proc == 'ww':
        pass

    elif proc == 'zz':
        pass

    else: 
        print("must enter a valid process (mumu, aa, ww, zz, vv). Quitting...")
        return 1
      
    # now write the lines to the file.. 
    with open(filename, 'w') as file:
        file.write('import model fermportal_QL\n')
        if (proc == 'vv' or proc == 'aa' or proc == 'zz' or proc == 'ww'):    
            file.write('set group_subprocesses false\n')
        if proc == 'vv': 
            file.write('define vxp = w+ z a\n')
            file.write('define vxm = w- z a\n')
        file.write(generate_line)
        file.write("output {}\n".format(run_name))

        # now we do the loops over masses (and other switches)
        if proc == 'mumu':
            # for mumu only, we have to loop over the PDF switches
            for switch1 in [0, 1]:
                for switch2 in [0, 1]:
                    for mass in mphi_values.tolist():
                        Wphiu = decaywidth_total_phiu(mass, 1e-6, 1e-7)
                        Wphid = decaywidth_total_phid(mass, 1e-6, 1e-7)
                    
                        file.write('launch {}\n'.format(run_name))
                        file.write('set run_tag mphi_{}_switch_{}_{}\n'.format(mass, switch1, switch2))
                        file.write(base_script)
                        file.write('set pdlabel lhapdf\n')
                        file.write('set lhaid {}\n'.format(resc_lhaid[pdf_resc]))
                        file.write('set lpp1 {}\n'.format(switch1))
                        file.write('set lpp2 {}\n'.format(switch2))
                        file.write('set fixed_fac_scale True\n')
                        if fact_scale == 'dyn':
                            file.write('set dsqrt_q2fact1 {}\n'.format(mass))
                            file.write('set dsqrt_q2fact2 {}\n'.format(mass))
                        else:
                            file.write('set dsqrt_q2fact1 5000.0\n')
                            file.write('set dsqrt_q2fact2 5000.0\n')
                        file.write('set msca {}\n'.format(mass))
                        file.write('set decay 9000006 {}\n'.format(Wphiu))
                        file.write('set decay 9000007 {}\n'.format(Wphid))

        else: 
            # for other processes, just loop over the masses
            for mass in mphi_values.tolist():
                Wphiu = decaywidth_total_phiu(mass, 1e-6, 1e-7)
                Wphid = decaywidth_total_phid(mass, 1e-6, 1e-7)
            
                file.write('launch {}\n'.format(run_name))
                file.write('set run_tag mphi_{}\n'.format(mass))
                file.write(base_script)

                file.write(eva_script)

                file.write('set msca {}\n'.format(mass))
                file.write('set decay 9000006 {}\n'.format(Wphiu))
                file.write('set decay 9000007 {}\n'.format(Wphid))
   
        # things to write to every file
        file.write('launch {} -i\n'.format(run_name))
        file.write('print_results --path=./cross_section_{}.txt --format=short\n'.format(run_name))


def execute_mg5():
    pass

#if __name__ == "__main__":
    
os.system('echo "Generating files ... "')

#create_script_file('mumu', 'L', 'R')
create_script_file('mumu', 'L', 'R', pdf_resc='1e-4', fact_scale='dyn')
create_script_file('mumu', 'R', 'L', pdf_resc='1e-4', fact_scale='dyn')

create_script_file('mumu', 'L', 'R', pdf_resc='1e-8', fact_scale='dyn')
create_script_file('mumu', 'R', 'L', pdf_resc='1e-8', fact_scale='dyn')

create_script_file('vv', None, None)
#create_script_file('ww', None, None)
#create_script_file('aa', None, None)


