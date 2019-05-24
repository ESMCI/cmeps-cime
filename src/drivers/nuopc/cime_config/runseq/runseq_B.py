#!/usr/bin/env python

import os, shutil, sys

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","..")
sys.path.append(os.path.join(_CIMEROOT, "scripts", "Tools"))

from standard_script_setup import *

logger = logging.getLogger(__name__)

def runseq(case, coupling_times):

    rundir    = case.get_value("RUNDIR")
    caseroot  = case.get_value("CASEROOT")
    comp_atm  = case.get_value("COMP_ATM")
    comp_ice  = case.get_value("COMP_ICE")
    comp_glc  = case.get_value("COMP_GLC")
    comp_lnd  = case.get_value("COMP_LND")
    comp_ocn  = case.get_value("COMP_OCN")
    comp_rof  = case.get_value("COMP_ROF")
    comp_wav  = case.get_value("COMP_WAV")
    budgets   = case.get_value("BUDGETS")

    glc_cpl_dt = coupling_times["glc_cpl_dt"]
    rof_cpl_dt = coupling_times["rof_cpl_dt"]
    ocn_cpl_dt = coupling_times["ocn_cpl_dt"]
    atm_cpl_dt = coupling_times["atm_cpl_dt"]

    outfile   = open(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), "w")

    # assume RASM_OPTION1

    outfile.write ("runSeq::                                 \n")
    if comp_glc == 'cism':
        outfile.write ("@" + str(glc_cpl_dt) + "             \n" ) # start of glc time loop
        outfile.write ("  MED med_phases_prep_glc            \n" )
        outfile.write ("  MED -> GLC :remapMethod=redist     \n" )

    if comp_rof == 'mosart' or comp_rof == 'rtm':
        outfile.write ("@" + str(rof_cpl_dt) + "             \n" ) # start of rof time loop
        outfile.write ("  MED med_phases_prep_rof_accum_avg  \n" )
        outfile.write ("  MED med_phases_prep_rof_avg        \n" )
        outfile.write ("  MED -> ROF :remapMethod=redist     \n" )

    outfile.write ("@" + str(ocn_cpl_dt) + "                 \n" ) # start of ocn time loop
    outfile.write ("  MED med_phases_prep_ocn_accum_avg      \n" )
    outfile.write ("  MED -> OCN :remapMethod=redist         \n" )
    outfile.write ("  OCN                                    \n" )

    outfile.write ("@" + str(atm_cpl_dt) + "                 \n" ) # start of atm time loop
    outfile.write ("  MED med_phases_prep_ocn_map            \n" )
    outfile.write ("  MED med_phases_aofluxes_run            \n" )
    outfile.write ("  MED med_phases_prep_ocn_merge          \n" )
    outfile.write ("  MED med_phases_prep_ocn_accum_fast     \n" )
    outfile.write ("  MED med_phases_ocnalb_run              \n" )
    if budgets:
        outfile.write ("  MED med_phases_diag_ocn            \n" )

    outfile.write ("  MED med_phases_prep_lnd                \n" )
    outfile.write ("  MED -> LND :remapMethod=redist         \n" )

    if budgets:
        outfile.write ("  MED med_phases_diag_lnd            \n" ) # budgets1 for cesm
        outfile.write ("  MED med_phases_diag_ice_med2ice    \n" )
        if comp_rof == 'mosart' or comp_rof == 'rtm':
            outfile.write ("  MED med_phases_diag_rof        \n" ) 

    outfile.write ("  MED med_phases_prep_ice                \n" )
    outfile.write ("  MED -> ICE :remapMethod=redist         \n" )

    if comp_wav == 'ww':
        outfile.write ("  MED med_phases_prep_wav            \n" )
        outfile.write ("  MED -> WAV :remapMethod=redist     \n" )

    if comp_wav != 'swav':
        outfile.write ("  WAV                                \n" )

    outfile.write ("  ICE                                    \n" )
    outfile.write ("  LND                                    \n" )
    outfile.write ("  LND -> MED :remapMethod=redist         \n" )
    if comp_rof == 'mosart' or comp_rof == 'rtm':
        outfile.write ("  MED med_phases_prep_rof_accum_fast \n" )

    outfile.write ("  ICE -> MED :remapMethod=redist         \n" )
    outfile.write ("  MED med_fraction_set                   \n" )

    outfile.write ("  MED med_phases_prep_atm                \n" )
    outfile.write ("  MED -> ATM :remapMethod=redist         \n" )
    outfile.write ("  ATM                                    \n" )
    outfile.write ("  ATM -> MED :remapMethod=redist         \n" )

    if comp_wav != 'swav':
        outfile.write ("  WAV -> MED :remapMethod=redist     \n" )

    if budgets:
        outfile.write ("  MED med_phases_diag_atm            \n" ) # budgets2 for cesm
        outfile.write ("  MED med_phases_diag_ice_ice2med    \n" )
        outfile.write ("  MED med_phases_diag_accum          \n" )
        outfile.write ("  MED med_phases_diag_print          \n" )

    outfile.write ("@                                        \n" ) # end of atm time loop

    outfile.write ("  OCN -> MED :remapMethod=redist         \n" )

    if comp_rof == 'mosart' or comp_rof == 'rtm':
        outfile.write ("@                                    \n" ) # end of ocn time loop
        outfile.write ("  ROF                                \n" )
        outfile.write ("  ROF -> MED :remapMethod=redist     \n" )

    if comp_glc == 'cism':
        outfile.write ("@                                    \n" ) # end of rof time loop
        outfile.write ("GLC                                  \n" )
        outfile.write ("                                     \n" )
        outfile.write ("GLC -> MED :remapMethod=redist       \n" )
  
    outfile.write ("  MED med_phases_history_write           \n" )
    outfile.write ("  MED med_phases_restart_write           \n" )
    outfile.write ("  MED med_phases_profile                 \n" )
    outfile.write ("@                                        \n" ) # end of glc or rof time loop
    outfile.write ("::                                       \n" )

    outfile.close()
    shutil.copy(os.path.join(caseroot, "CaseDocs", "nuopc.runseq"), rundir)
