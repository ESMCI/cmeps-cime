module med_phases_diag_mod

contains

  subroutine cime_run_calc_budgets1()

    !----------------------------------------------------------
    ! Budget with old fractions
    !----------------------------------------------------------

    ! WJS (2-17-11): I am just using the first instance for the budgets because we
    ! don't expect budgets to be conserved for our case (I case). Also note that we
    ! don't expect budgets to be conserved for the interactive ensemble use case either.
    ! tcraig (aug 2012): put this after rof->cpl so the budget sees the new r2x_rx.
    ! it will also use the current r2x_ox here which is the value from the last timestep
    ! consistent with the ocean coupling

    call t_drvstartf ('CPL:BUDGET1',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
    if (is_local%wrap%comp_present(complnd)) then 

       call med_diag_lnd(is_local%wrap%FBImp(complnd,complnd), is_local%wrap%FBExp(complnd), &
            is_local%wrap%FBfrac(complnd), do_l2x=.true., do_x2l=.true.)
    endif
    if (is_local%wrap%comp_present(comprof)) then 
       call med_diag_rof(rof, fractions_rx, infodata)
    endif
    if (is_local%wrap%comp_present(compice)) then 
       call med_diag_ice(ice, fractions_ix, infodata, do_x2i=.true.)
    endif
    call t_drvstopf  ('CPL:BUDGET1',cplrun=.true.,budget=.true.)

  end subroutine cime_run_calc_budgets1

!----------------------------------------------------------------------------------

  subroutine cime_run_calc_budgets2()

    !----------------------------------------------------------
    ! Budget with new fractions
    !----------------------------------------------------------

       call t_drvstartf ('CPL:BUDGET2',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
       if (atm_present) then
          call med_diag_atm(atm, fractions_ax, infodata, do_a2x=.true., do_x2a=.true.)
       endif
       if (ice_present) then
          call med_diag_ice(ice, fractions_ix, infodata, do_i2x=.true.)
       endif
       call t_drvstopf  ('CPL:BUDGET2',cplrun=.true.,budget=.true.)

       call t_drvstartf ('CPL:BUDGET3',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
       call med_diag_accum()
       call t_drvstopf  ('CPL:BUDGET3',cplrun=.true.,budget=.true.)

       call t_drvstartf ('CPL:BUDGETF',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
       if (.not. dead_comps) then
          call med_diag_print(EClock_d,stop_alarm,budget_inst, &
               budget_daily, budget_month, budget_ann, budget_ltann, &
               budget_ltend, infodata)
       endif
       call med_diag_zero(EClock=EClock_d)

       call t_drvstopf  ('CPL:BUDGETF',cplrun=.true.,budget=.true.)
    end if
  end subroutine cime_run_calc_budgets2

!----------------------------------------------------------------------------------

  subroutine cime_run_calc_budgets3()

    !----------------------------------------------------------
    ! ocn budget (rasm_option2)
    !----------------------------------------------------------

    if (iamin_CPLID) then
       call cime_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
       call t_drvstartf ('CPL:BUDGET0',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
       xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
       call med_diag_ocn(ocn, xao_ox(1), fractions_ox, infodata, &
            do_o2x=.true., do_x2o=.true., do_xao=.true.)
       call t_drvstopf ('CPL:BUDGET0',cplrun=.true.,budget=.true.)
    end if
  end subroutine cime_run_calc_budgets3

end module med_phases_diag_mod
