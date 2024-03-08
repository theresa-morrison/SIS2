!> Provides a common interface for jointly stepping SIS2 and MOM6, and will
!! evolve as a platform for tightly integrating the ocean and sea ice models.
module SIS_dyn_setup

! This file is a part of SIS2. See LICENSE.md for the license.

!-----------------------------------------------------------------------
!
! This module provides a common interface for jointly stepping SIS2 and MOM6, and
! will evolve as a platform for tightly integrating the ocean and sea ice models.

use MOM_error_handler, only : SIS_error=>MOM_error, FATAL, WARNING, SIS_mesg=>MOM_mesg
use MOM_time_manager,   only : time_type, time_type_to_real 
use MOM_time_manager,   only : operator(+), operator(-), operator(>) 

use ice_model_mod,      only : ice_data_type 

use SIS_types,         only : IST_chksum, IST_bounds_check
use ice_grid,          only : ice_grid_type
use SIS_hor_grid,      only : SIS_hor_grid_type 
use SIS_types,         only : ice_state_type
use SIS_types,         only : fast_ice_avg_type
use MOM_unit_scaling,  only : unit_scale_type
use SIS_types,         only : ocean_sfc_state_type
use SIS_types,         only : ice_ocean_flux_type
use SIS_dyn_trans,     only : dyn_trans_CS
!use ice_bergs,         only : iceberg
!use SIS_tracer_flow_control, only : SIS_tracer_flow_control_CS
use SIS_open_boundary, only : ice_OBC_type
use SIS_dyn_trans,     only : dyn_State_2d
use MOM_forcing_type,  only : SIS_C_EVP_state
use MOM_SIS_C_dyn_CS_type, only : SIS_C_dyn_CS

use SIS_dyn_trans, only : convert_IST_to_simple_state, increase_max_tracer_step_memory  
use SIS_dyn_trans, only : set_wind_stresses_C

use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,     only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP, CLOCK_ROUTINE
use MOM_domains,       only : pass_var, pass_vector, AGRID, BGRID_NE, CGRID_NE
use MOM_time_manager,  only : operator(+), operator(-)
use MOM_time_manager,  only : operator(>), operator(*), operator(/), operator(/=)
use MOM_unit_scaling,  only : unit_scale_type

use SIS_debugging,     only : chksum, Bchksum, hchksum
use SIS_debugging,     only : hchksum_pair, Bchksum_pair, uvchksum
use SIS_diag_mediator, only : enable_SIS_averaging, disable_SIS_averaging
use SIS_diag_mediator, only : query_SIS_averaging_enabled, SIS_diag_ctrl
use SIS_dyn_cgrid,     only : basal_stress_coeff_C
use SIS_dyn_cgrid,     only : basal_stress_coeff_itd
use SIS_dyn_cgrid,     only : limit_stresses 
 
use ice_type_mod,       only : Ice_public_type_chksum, Ice_public_type_bounds_check
use MOM_time_manager,  only : time_type, time_type_to_real, real_to_time

use MOM_domains,       only : CORNER 


#include <SIS2_memory.h>

implicit none ; private

public :: setup_SIS_dynamics 

integer :: iceClock
integer :: ice_clock_slow, ice_clock_fast, ice_clock_exchange

contains

!> This subroutine is used to call SIS dynamics directly from MOM so it needs to public and visible to MOM
! update_ice_dynamics_trans(Ice, time_step=dyn_time_step, start_cycle=(ns==1), end_cycle=(ns==nstep), cycle_length=dt_coupling)
subroutine setup_SIS_dynamics(Ice, EVPT, time_step, cycle_length)
  type(ice_data_type),       intent(inout) :: Ice !< The publicly visible ice data type.
  type(SIS_C_EVP_state),     intent(inout) :: EVPT !< 
  type(time_type), optional, intent(in)    :: time_step !< The amount of time to cover in this update.
  real           , optional, intent(in)    :: cycle_length
 
  ! These pointers are used to simplify the code below.
  type(ice_grid_type),              pointer :: sIG => NULL()
  type(SIS_hor_grid_type),          pointer :: sG => NULL()
  type(ice_state_type),             pointer :: sIST => NULL()
  type(fast_ice_avg_type),          pointer :: FIA => NULL()
  type(unit_scale_type),            pointer :: US => NULL()
  type(ocean_sfc_state_type),       pointer :: OSS => NULL() !< A structure containing the arrays that describe
  type(ice_ocean_flux_type),        pointer :: IOF => NULL() !< A structure containing fluxes from the ice to
  type(dyn_trans_CS),               pointer :: CS =>  NULL() !< The control structure for the SIS_dyn_trans module
  ! type(icebergs),                   pointer :: icebergs_CS => NULL() !< A control structure for the iceberg model.
  ! type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp => NULL()  !< The structure for controlling calls to
  type(ice_OBC_type),               pointer :: OBC => NULL() !< Open boundary structure.
  type(dyn_state_2d),               pointer :: DS2d => NULL() !< A simplified 2-d description of the ice state

  real :: dt_slow  ! The time step over which to advance the model [T ~> s].
  logical :: do_multi_trans 

  ! Local variables from SIS_multi_dyn_trans
  real :: dt_adv_cycle ! The length of the advective cycle timestep [T ~> s].
  real :: dt_diags     ! The length of time over which the diagnostics are valid [T ~> s].
  type(time_type) :: Time_cycle_start ! The model's time at the start of an advective cycle.
  integer :: nadv_cycle, nac ! The number of tracer advective cycles within this call.
  ! Local variables SIS_merged_dyn_cont
  real, dimension(SZI_(Ice%sCS%G),SZJ_(Ice%sCS%G))   :: &
    ice_free, &         ! The fractional open water [nondim], between 0 & 1.
    rdg_rate            ! A ridging rate [T-1 ~> s-1], this will be calculated from the strain rates
                        ! in the dynamics.
  real, dimension(SZIB_(Ice%sCS%G),SZJB_(Ice%sCS%G)) :: &
    WindStr_x_B, &      ! Zonal (_x_) and meridional (_y_) wind stresses
    WindStr_y_B, &      ! averaged over the ice categories on a B-grid [R Z L T-2 ~> Pa].
    WindStr_x_ocn_B, &  ! Zonal wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    WindStr_y_ocn_B, &  ! Meridional wind stress on the ice-free ocean on a B-grid [R Z L T-2 ~> Pa].
    str_x_ice_ocn_B, &  ! Zonal ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
    str_y_ice_ocn_B     ! Meridional ice-ocean stress on a B-grid [R Z L T-2 ~> Pa].
  real, dimension(SZIB_(Ice%sCS%G),SZJ_(Ice%sCS%G))  :: &
    WindStr_x_Cu, &   ! Zonal wind stress averaged over the ice categories on C-grid u-points [R Z L T-2 ~> Pa].
    WindStr_x_ocn_Cu, & ! Zonal wind stress on the ice-free ocean on C-grid u-points [R Z L T-2 ~> Pa].
    str_x_ice_ocn_Cu   ! Zonal ice-ocean stress on C-grid u-points [R Z L T-2 ~> Pa].
  real, dimension(SZI_(Ice%sCS%G),SZJB_(Ice%sCS%G))  :: &
    WindStr_y_Cv, &   ! Meridional wind stress averaged over the ice categories on C-grid v-points [R Z L T-2 ~> Pa].
    WindStr_y_ocn_Cv, & ! Meridional wind stress on the ice-free ocean on C-grid v-points [R Z L T-2 ~> Pa].
    str_y_ice_ocn_Cv  ! Meridional ice-ocean stress on C-grid v-points [R Z L T-2 ~> Pa].
  
  real :: ps_vel   ! The fractional thickness category coverage at a velocity point.
  real :: wt_new, wt_prev ! Weights in an average.
  real :: dt_slow_dyn  ! The slow dynamics timestep [T ~> s].
  real :: dt_slow_dyn_sec ! The slow dynamics timestep [s].
  real :: dt_adv       ! The advective subcycle timestep [T ~> s].
  logical :: continuing_call ! If true, there are more in the series of advective updates
                             ! after this call.
  integer :: ndyn_steps, nds ! The number of dynamic steps in this call.
  integer :: i, j, k, n, isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed !, IsdB, IedB, JsdB, JedB

  logical :: cycle_start, cycle_end, end_of_cycle
  
  ! ..... SIS_dynamics_trans........................................................
  if (.not.associated(Ice%sCS)) call SIS_error(FATAL, &
      "The pointer to Ice%sCS must be associated in update_ice_dynamics_trans.")

  sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G ; FIA => Ice%sCS%FIA ; US => Ice%sCS%US
  OSS => Ice%sCS%OSS ; IOF => Ice%sCS%IOF ;  ! icebergs_CS => Ice%icebergs ; ! tracer_CSp => Ice%sCS%SIS_tracer_flow_CSp
  OBC => Ice%OBC ;  CS => Ice%sCS%dyn_trans_CSp ; DS2d => Ice%sCS%dyn_trans_CSp%DS2d
  
  if (.not.CS%merged_cont) call SIS_error(FATAL, &
          "SIS_multi_dyn_trans should not be called unless MERGED_CONTINUITY=True.")

  dt_slow = US%s_to_T*time_type_to_real(Ice%sCS%Time_step_slow)
  if (present(time_step)) dt_slow = US%s_to_T*time_type_to_real(time_step)

  call cpu_clock_begin(iceClock) ; call cpu_clock_begin(ice_clock_slow)

  ! Do halo updates on the forcing fields, as necessary.  This must occur before
  ! the call to SIS_dynamics_trans, because update_icebergs does its own halo
  ! updates, and slow_thermodynamics only works on the computational domain.
  call pass_vector(FIA%WindStr_x, FIA%WindStr_y, sG%Domain, stagger=AGRID, complete=.false.)
  call pass_vector(FIA%WindStr_ocn_x, FIA%WindStr_ocn_y, sG%Domain, stagger=AGRID)
  call pass_var(FIA%ice_cover, sG%Domain, complete=.false.)
  call pass_var(FIA%ice_free,  sG%Domain, complete=.true.)
  if (sIST%valid_IST) then
    call pass_var(sIST%part_size, sG%Domain)
    call pass_var(sIST%mH_ice, sG%Domain, complete=.false.)
    call pass_var(sIST%mH_pond, sG%Domain, complete=.false.)
    call pass_var(sIST%mH_snow, sG%Domain, complete=.true.)
  endif

  if (Ice%sCS%debug) then
    call Ice_public_type_chksum("Before SIS_dynamics_trans", Ice, check_slow=.true.)
  endif

  ! Now ............................................................................................
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !> SIS_multi_dyn_trans makes the calls to do ice dynamics and mass and tracer transport as
  !! appropriate for a dynamic and advective update cycle with multiple calls.
  !subroutine SIS_multi_dyn_trans(IST, OSS, FIA, IOF, dt_slow, CS, icebergs_CS, G, US, IG, tracer_CSp, &
  !                               OBC, start_cycle, end_cycle, cycle_length)
  IOF%stress_count = 0

  dt_diags = dt_slow ; if (present(cycle_length)) dt_diags = US%s_to_T*cycle_length
  nadv_cycle = 1; dt_adv_cycle = dt_slow / real(nadv_cycle)

  ! Convert the category-resolved ice state into the simplified 2-d ice state.
  ! This should be called after a thermodynamic step or if ice_transport was called.
  call convert_IST_to_simple_state(sIST, CS%DS2d, CS%CAS, sG, US, sIG, CS)

  ! Update the category-merged dynamics and use the merged continuity equation.
  ! This could be called as many times as necessary.
  Time_cycle_start = CS%Time - real_to_time((nadv_cycle-(nac-1))*US%T_to_s*dt_adv_cycle)
  end_of_cycle = (nac < nadv_cycle) .or. cycle_end
  ! Update the category-merged dynamics and use the merged continuity equation.
  ! call SIS_merged_dyn_cont(OSS, FIA, IOF, CS%DS2d, IST, dt_adv_cycle, Time_cycle_start, &
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  !> Update the category-merged ice state and call the merged continuity update.
  !subroutine SIS_merged_dyn_cont(OSS, FIA, IOF, DS2d, IST, dt_cycle, Time_start, G, US, IG, CS, OBC, end_call)
  ! This subroutine updates the 2-d sea-ice dynamics.
  !   Variables updated here: DS2d%ice_cover, DS2d%[uv]_ice_[BC], DS2d%mca_step, DS2d%mi_sum,
  !       CS%[uv]h_step, DS2d%nts, CS%SIS_[BC]_dyn_CSp,  IOF (stresses)
    
  isc = sG%isc ; iec = sG%iec ; jsc = sG%jsc ; jec = sG%jec
  isd = sG%isd ; ied = sG%ied ; jsd = sG%jsd ; jed = sG%jed

  ndyn_steps = 1
  dt_slow_dyn = dt_adv_cycle / ndyn_steps
  dt_slow_dyn_sec = US%T_to_s*dt_slow_dyn
  dt_adv = dt_slow_dyn / real(CS%adv_substeps)
  if (ndyn_steps*CS%adv_substeps > DS2d%max_nts) &
    call increase_max_tracer_step_memory(DS2d, sG, ndyn_steps*CS%adv_substeps)
  continuing_call = .false.  ! ; if (present(end_call)) continuing_call = .not.end_call

  !call cpu_clock_begin(iceClock4)
  call enable_SIS_averaging(dt_slow_dyn_sec, Time_cycle_start + real_to_time(nds*dt_slow_dyn_sec), CS%diag)
  do j=jsd,jed ; do i=isd,ied ; ice_free(i,j) = max(1.0 - DS2d%ice_cover(i,j), 0.0) ; enddo ; enddo

  ! In the dynamics code, only the ice velocities are changed, and the ice-ocean
  ! stresses are calculated.  The gravity wave dynamics (i.e. the continuity
  ! equation) are not included in the dynamics (yet).  All of the thickness categories
  ! are merged together.

  ! Correct the wind stresses for changes in the fractional ice-coverage and set
  ! the wind stresses on the ice and the open ocean for a C-grid staggering.
  ! This block of code must be executed if ice_cover and ice_free or the various wind
  ! stresses were updated.
  call set_wind_stresses_C(FIA, DS2d%ice_cover, ice_free, WindStr_x_Cu, WindStr_y_Cv, &
                           WindStr_x_ocn_Cu, WindStr_y_ocn_Cv, sG, US, CS%complete_ice_cover, OBC)

  if (CS%lemieux_landfast) then
    call basal_stress_coeff_C(sG, DS2d%mi_sum, DS2d%ice_cover, OSS%sea_lev, CS%SIS_C_dyn_CSp)
  elseif (CS%itd_landfast) then
    call basal_stress_coeff_itd(sG, sIG, sIST, OSS%sea_lev, CS%SIS_C_dyn_CSp)
  endif

  if (CS%debug) then
    call uvchksum("Before SIS_C_dynamics [uv]_ice_C", DS2d%u_ice_C, DS2d%v_ice_C, sG, scale=US%L_T_to_m_s)
    call hchksum(ice_free, "ice_free before SIS_C_dynamics", sG%HI)
    call hchksum(DS2d%mca_step(:,:,DS2d%nts), "misp_sum before SIS_C_dynamics", sG%HI, scale=US%RZ_to_kg_m2)
    call hchksum(DS2d%mi_sum, "mi_sum before SIS_C_dynamics", sG%HI, scale=US%RZ_to_kg_m2)
    call hchksum(OSS%sea_lev, "sea_lev before SIS_C_dynamics", sG%HI, haloshift=1, scale=US%Z_to_m)
    call hchksum(DS2d%ice_cover, "ice_cover before SIS_C_dynamics", sG%HI, haloshift=1)
    call uvchksum("[uv]_ocn before SIS_C_dynamics", OSS%u_ocn_C, OSS%v_ocn_C, sG, halos=1, scale=US%L_T_to_m_s)
    call uvchksum("WindStr_[xy] before SIS_C_dynamics", WindStr_x_Cu, WindStr_y_Cv, sG, &
                  halos=1, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
    !call hchksum_pair("WindStr_[xy]_A before SIS_C_dynamics", WindStr_x_A, WindStr_y_A, G, halos=1)
  endif

  !call cpu_clock_begin(iceClocka)
  !### Ridging needs to be added with C-grid dynamics.
  if (CS%do_ridging) rdg_rate(:,:) = 0.0
  call SIS_C_dyn_to_EVP(DS2d%ice_cover, DS2d%mca_step(:,:,DS2d%nts), DS2d%mi_sum, &
                        DS2d%u_ice_C, DS2d%v_ice_C, &
                        OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
                        str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, sG, US, CS%SIS_C_dyn_CSp, EVPT)

  ! back to SIS_dynamics_trans ..............................................................................
  ! Complete the category-resolved mass and tracer transport and update the ice state type.
  ! This should go - at least after the ocean dynamics ...................................
  ! call complete_IST_transport(CS%DS2d, CS%CAS, IST, dt_adv_cycle, G, US, IG, CS)

  !  if (CS%column_check .and. (DS2d%nts==0)) &
  !    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
  !                              message="      Post_transport")! , check_column=.true.)

  ! Finalized the streses for use by the ocean.
  ! call finish_ocean_top_stresses(IOF, G)

  ! Do diagnostics and update some information for the atmosphere.
  ! call ice_state_cleanup(IST, OSS, IOF, dt_slow, G, US, IG, CS, tracer_CSp)
 
  ! back to updat_ice_dynamics.......................................................................
  ! Set up the stresses and surface pressure in the externally visible structure Ice.
  ! if (sIST%valid_IST) call ice_mass_from_IST(sIST, Ice%sCS%IOF, sG, sIG)

  ! if (Ice%sCS%debug) then
  !  call IOF_chksum("Before set_ocean_top_dyn_fluxes", Ice%sCS%IOF, sG, US, mech_fluxes=.true.)
  ! endif
  ! call set_ocean_top_dyn_fluxes(Ice, Ice%sCS%IOF, FIA, sG, US, Ice%sCS)

end subroutine setup_SIS_dynamics

!call SIS_C_dynamics(DS2d%ice_cover, DS2d%mca_step(:,:,DS2d%nts), DS2d%mi_sum, &
!                    DS2d%u_ice_C, DS2d%v_ice_C, &
!                    OSS%u_ocn_C, OSS%v_ocn_C, WindStr_x_Cu, WindStr_y_Cv, OSS%sea_lev, &
!                    str_x_ice_ocn_Cu, str_y_ice_ocn_Cv, dt_slow_dyn, G, US, CS%SIS_C_dyn_CSp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!> SIS_C_dynamics takes a single dynamics timestep with EVP subcycles
subroutine SIS_C_dyn_to_EVP(ci, mis, mice, ui, vi, uo, vo, fxat, fyat, &
                         sea_lev, fxoc, fyoc, dt_slow, G, US, CS, EVPT)
     type(SIS_hor_grid_type),           intent(inout) :: G   !< The horizontal grid type
     real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: ci  !< Sea ice concentration [nondim]
     real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: mis   !< Mass per unit ocean area of sea ice,
                                                               !! snow and melt pond water [R Z ~> kg m-2]
     real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
     real, dimension(SZIB_(G),SZJ_(G)), intent(inout) :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
     real, dimension(SZI_(G),SZJB_(G)), intent(inout) :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]
     real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: uo    !< Zonal ocean velocity [L T-1 ~> m s-1]
     real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: vo    !< Meridional ocean velocity [L T-1 ~> m s-1]
     real, dimension(SZIB_(G),SZJ_(G)), intent(in   ) :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
     real, dimension(SZI_(G),SZJB_(G)), intent(in   ) :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]
     real, dimension(SZI_(G),SZJ_(G)),  intent(in   ) :: sea_lev !< The height of the sea level, including

     real, dimension(SZIB_(G),SZJ_(G)), intent(  out) :: fxoc  !< Zonal ice stress on ocean [R Z L T-2 ~> Pa]
     real, dimension(SZI_(G),SZJB_(G)), intent(  out) :: fyoc  !< Meridional ice stress on ocean [R Z L T-2 ~> Pa]
     real,                              intent(in   ) :: dt_slow !< The amount of time over which the ice
                                                               !! dynamics are to be advanced [T ~> s].
     type(unit_scale_type),             intent(in)    :: US    !< A structure with unit conversion factors
     type(SIS_C_dyn_CS),                pointer       :: CS    !< The control structure for this module
!     logical,         optional, intent(in)    :: first_call
!     logical,         optional, intent(in)    :: second_call
     type(SIS_C_EVP_state), &
                      optional, intent(inout) :: EVPT !! type containing the variables needed to do the EVP TJM
   
     ! Local variables
     real, dimension(SZI_(G),SZJ_(G)) :: &
       sh_Dt, &    ! sh_Dt is the horizontal tension (du/dx - dv/dy) including
                   ! all metric terms [T-1 ~> s-1].
       sh_Dd       ! sh_Dd is the flow divergence (du/dx + dv/dy) including all
                   ! metric terms [T-1 ~> s-1].
     real, dimension(SZIB_(G),SZJB_(G)) :: &
       sh_Ds       ! sh_Ds is the horizontal shearing strain (du/dy + dv/dx)
                   ! including all metric terms [T-1 ~> s-1].
   
   
     real, dimension(SZI_(G),SZJ_(G)) :: &
       pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
       ci_proj, &  ! The projected ice concentration [nondim].
       zeta, &     ! The ice bulk viscosity [R Z L2 T-1 ~> Pa m s] (i.e., [N s m-1]).
       del_sh, &   ! The magnitude of the shear rates [T-1 ~> s-1].
       diag_val, & ! A temporary diagnostic array.
       del_sh_min_pr, &  ! When multiplied by pres_mice, this gives the minimum
                   ! value of del_sh that is used in the calculation of zeta [T-1 ~> s-1].
                   ! This is set based on considerations of numerical stability,
                   ! and varies with the grid spacing.
       dx2T, dy2T, &   ! dx^2 or dy^2 at T points [L2 ~> m2].
       dx_dyT, dy_dxT, &  ! dx/dy or dy_dx at T points [nondim].
       siu, siv, sispeed  ! diagnostics on T points [L T-1 ~> m s-1].
   
     real, dimension(SZIB_(G),SZJ_(G)) :: &
       fxic, &   ! Zonal force due to internal stresses [R Z L T-2 ~> Pa].
       fxic_d, & ! Zonal force due to divergence internal stress [R Z L T-2 ~> Pa].
       fxic_t, & ! Zonal force due to tension internal stress [R Z L T-2 ~> Pa].
       fxic_s, & ! Zonal force due to shearing internal stress [R Z L T-2 ~> Pa].
       fxlf, &   ! Zonal landfast ice stress [R Z L T-2 ~> Pa]
       ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
       ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
       Cor_u, & ! Zonal Coriolis acceleration [L T-2 ~> m s-2].
       PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
       diag_val_u, & ! A temporary diagnostic array.
       u_tmp, & ! A temporary copy of the old values of ui [L T-1 ~> m s-1].
       u_IC, &  ! The initial zonal ice velocities [L T-1 ~> m s-1].
       mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
       f2dt_u, &! The squared effective Coriolis parameter at u-points times a
                ! time step [T-1 ~> s-1].
       I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].
     real, dimension(SZI_(G),SZJB_(G)) :: &
       fyic, &   ! Meridional force due to internal stresses [R Z L T-2 ~> Pa].
       fyic_d, & ! Meridional force due to divergence internal stress [R Z L T-2 ~> Pa].
       fyic_t, & ! Meridional force due to tension internal stress [R Z L T-2 ~> Pa].
       fyic_s, & ! Meridional force due to shearing internal stress [R Z L T-2 ~> Pa].
       fylf, &   ! Meridional landfast ice stress [R Z L T-2 ~> Pa]
       vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
       vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
       Cor_v, &  ! Meridional Coriolis acceleration [L T-2 ~> m s-2].
       PFv, &   ! Meridional hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
       diag_val_v, & ! A temporary diagnostic array.
       v_IC, &  ! The initial meridional ice velocities [L T-1 ~> m s-1].
       mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
       f2dt_v, &! The squared effective Coriolis parameter at v-points times a
                ! time step [T-1 ~> s-1].
       I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].
   
     real, dimension(SZIB_(G),SZJB_(G)) :: &
       mi_ratio_A_q, & ! A ratio of the masses interpolated to the faces around a
                ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
                ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].
       q, &     ! A potential-vorticity-like field for the ice, the Coriolis parameter
                ! divided by a spatially averaged mass per unit area [T-1 R-1 Z-1 ~> s-1 m2 kg-1].
       dx2B, dy2B, &   ! dx^2 or dy^2 at B points [L2 ~> m2].
       dx_dyB, dy_dxB  ! dx/dy or dy_dx at B points [nondim].
     real, dimension(SZIB_(G),SZJ_(G)) :: &
       azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
       czon, dzon, & ! are applied to the neighboring values of vi & ui,
       amer, bmer, & ! respectively to get the barotropic inertial rotation,
       cmer, dmer    ! in units of [T-1 ~> s-1].  azon and amer couple the same pair of
                     ! velocities, but with the influence going in opposite
                     ! directions.
   
     real :: Cor       ! A Coriolis acceleration [L T-2 ~> m s-2].
     real :: fxic_now  ! Zonal ice internal stress convergence [R Z L T-2 ~> Pa].
     real :: fyic_now  ! Meridional ice internal stress convergence [R Z L T-2 ~> Pa].
     real :: drag_u, drag_v ! Drag rates with the ocean at u & v points [R Z T-1 ~> kg m-2 s-1].
     real :: drag_LFu  ! Drag rates to the land for landfast ice at u points [R Z T-1 ~> kg m-2 s-1].
     real :: drag_LFv  ! Drag rates to the land for landfast ice at v points [R Z T-1 ~> kg m-2 s-1].
     real :: drag_max  ! A maximum drag rate allowed in the ocean [R Z T-1 ~> kg m-2 s-1].
     real :: tot_area  ! The sum of the area of the four neighboring cells [L2 ~> m2].
     real :: dxharm    ! The harmonic mean of the x- and y- grid spacings [L ~> m].
     real :: muq2, mvq2  ! The product of the u- and v-face masses per unit cell
                         ! area surrounding a vorticity point [R2 Z2 ~> kg2 m-4].
     real :: muq, mvq    ! The u- and v-face masses per unit cell area extrapolated
                         ! to a vorticity point on the coast [R Z ~> kg m-2].
     real :: min_rescale ! The smallest of the 4 surrounding values of rescale [nondim].
     real :: I_1pdt_T    ! 1.0 / (1.0 + dt_2Tdamp) [nondim].
     real :: I_1pE2dt_T  ! 1.0 / (1.0 + EC^2 * dt_2Tdamp) [nondim].
   
     real :: v2_at_u     ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
     real :: u2_at_v     ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
     real :: v2_at_u_min ! The squared v-velocity interpolated to u points [L2 T-2 ~> m2 s-2].
     real :: u2_at_v_min ! The squared u-velocity interpolated to v points [L2 T-2 ~> m2 s-2].
     real :: uio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
     real :: vio_init    ! Ice-ocean velocity differences [L T-1 ~> m s-1]
     real :: m_uio_explicit ! Ice-ocean x-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
     real :: m_vio_explicit ! Ice-ocean y-velocity differences times the ice mass [R Z L T-1 ~> kg m-1 s-1]
     real :: uio_pred    ! Ice-ocean x-velocity differences [L T-1 ~> m s-1]
     real :: vio_pred    ! Ice-ocean y-velocity differences [L T-1 ~> m s-1]
     real :: I_cdRhoDt   ! The inverse of the product of the drag coefficient, ocean density and
                         ! timestep [L Z-1 R-1 T-1 ~> m3 kg-1 s-1].
     real :: cdRho       ! The ice density times the drag coefficient and rescaling factors [R Z L-1 ~> kg m-3]
     real :: b_vel0      ! The initial difference between the velocity magnitude
                         ! and the absolute value of the u- or v- component, plus
                         ! the ice thickness divided by the time step and the drag
                         ! coefficient [L T-1 ~> m s-1].
     real :: uio_C   ! A u-velocity difference between the ocean and ice [L T-1 ~> m s-1].
     real :: vio_C   ! A v-velocity difference between the ocean and ice [L T-1 ~> m s-1].
   
     real :: Tdamp   ! The damping timescale of the stress tensor components
                     ! toward their equilibrium solution due to the elastic terms [T ~> s].
     real :: dt      ! The short timestep associated with the EVP dynamics [T ~> s].
     real :: dt_2Tdamp ! The ratio of the timestep to the elastic damping timescale [nondim].
     real :: dt_cumulative ! The elapsed time within this call to EVP dynamics [T ~> s].
     integer :: EVP_steps ! The number of EVP sub-steps that will actually be taken.
     real :: I_sub_steps  ! The number inverse of the number of EVP time steps per
                     ! slow time step.
     real :: EC2     ! EC^2, where EC is the yield curve axis ratio.
     real :: I_EC2   ! 1/EC^2, where EC is the yield curve axis ratio.
     real :: I_EC    ! 1/EC, where EC is the yield curve axis ratio.
     real :: I_2EC   ! 1/(2*EC), where EC is the yield curve axis ratio.
     real, parameter :: H_subroundoff = 1e-30 ! A negligible ice thickness [m].
     real :: m_neglect  ! A tiny mass per unit area [R Z ~> kg m-2].
     real :: m_neglect2 ! A tiny mass per unit area squared [R2 Z2 ~> kg2 m-4].
     real :: m_neglect4 ! A tiny mass per unit area to the 4th power [R4 Z4 ~> kg4 m-8].
     real :: sum_area   ! The sum of ocean areas around a vorticity point [L2 ~> m2].
   
     type(time_type) :: &
       time_it_start, &  ! The starting time of the iterative steps.
       time_step_end, &  ! The end time of an iterative step.
       time_end_in       ! The end time for diagnostics when this routine started.
     real :: time_int_in ! The diagnostics' time interval when this routine started.
     logical :: do_hifreq_output  ! If true, output occurs every iterative step.
     logical :: do_trunc_its  ! If true, overly large velocities in the iterations are truncated.
     integer :: halo_sh_Ds  ! The halo size that can be used in calculating sh_Ds.
     integer :: i, j, isc, iec, jsc, jec, n
     isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
   
     if (.not.associated(CS)) call SIS_error(FATAL, &
            "SIS_C_dynamics: Module must be initialized before it is used.")
   
     if ((isc - G%isdB < 2) .or. (jsc - G%jsdB < 2)) call SIS_error(FATAL, &
            "SIS_C_dynamics is written to require a 2-point halo or 1-point and symmetric memory.")
   
     halo_sh_Ds = min(isc-G%isd, jsc-G%jsd, 2)
   
     ! Zero these arrays to accumulate sums.
     fxoc(:,:) = 0.0 ; fyoc(:,:) = 0.0
     fxlf(:,:) = 0.0 ; fylf(:,:) = 0.0
     fxic(:,:) = 0.0 ; fyic(:,:) = 0.0
     Cor_u(:,:) = 0.0 ; Cor_v(:,:) = 0.0
     fxic_d(:,:) = 0.0 ; fyic_d(:,:) = 0.0 ; fxic_t(:,:) = 0.0 ; fyic_t(:,:) = 0.0
     fxic_s(:,:) = 0.0 ; fyic_s(:,:) = 0.0
   
     if ((CS%evp_sub_steps<=0) .and. (CS%dt_Rheo<=0.0)) return
   
     if (CS%FirstCall) then
       !   These halo updates are only needed if the str_... arrays have just
       ! been read from a restart file.  Otherwise the halos contain valid data.
       call pass_var(CS%str_d, G%Domain) ; call pass_var(CS%str_t, G%Domain)
       call pass_var(CS%str_s, G%Domain, position=CORNER)
       CS%FirstCall = .false.
     endif
   
     if (CS%dt_Rheo > 0.0) then
       EVP_steps = max(CEILING(dt_slow/CS%dt_Rheo - 0.0001), 1)
     else
       EVP_steps = CS%evp_sub_steps
     endif
     dt = dt_slow/EVP_steps
   
     drag_max = CS%Rho_ocean * CS%min_ocn_inertial_h / dt_slow
     I_cdRhoDt = 1.0 / (CS%cdw * US%L_to_Z*CS%Rho_ocean * dt)
     do_trunc_its = (CS%CFL_check_its .and. (CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0))
   
     EC2 = CS%EC**2
     I_EC = 0.0 ; if (CS%EC > 0.0) I_EC = 1.0 / CS%EC
     I_2EC = 0.0 ; if (CS%EC > 0.0) I_2EC = 0.5 / CS%EC
     I_EC2 = 0.0 ; if (EC2 > 0.0) I_EC2 = 1.0 / EC2
   
     do_hifreq_output = .false.
     if ((CS%id_ui_hifreq > 0) .or. (CS%id_vi_hifreq > 0) .or. &
         (CS%id_str_d_hifreq > 0) .or. (CS%id_str_t_hifreq > 0) .or. &
         (CS%id_str_s_hifreq > 0) .or. (CS%id_sh_d_hifreq > 0) .or. &
         (CS%id_sh_t_hifreq > 0) .or. (CS%id_sh_s_hifreq > 0) .or. &
         (CS%id_ci_hifreq > 0) .or. (CS%id_stren_hifreq > 0)) then
       do_hifreq_output = query_SIS_averaging_enabled(CS%diag, time_int_in, time_end_in)
       if (do_hifreq_output) &
         time_it_start = time_end_in - real_to_time(US%T_to_s*dt_slow)
     endif
   
     Tdamp = CS%Tdamp
     if (CS%Tdamp == 0.0) then
       ! Hunke (2001) chooses a specified multiple (0.36) of dt_slow for Tdamp, and shows that
       ! stability requires Tdamp > 2*dt.  Here 0.2 is used instead for greater stability.
       Tdamp = max(0.2*dt_slow, 3.0*dt)
     elseif (CS%Tdamp < 0.0) then
       Tdamp = max(-CS%Tdamp*dt_slow, 3.0*dt)
     endif
     dt_2Tdamp = dt / (2.0 * Tdamp)
   
     ui_min_trunc(:,:) = 0.0 ; ui_max_trunc(:,:) = 0.0
     vi_min_trunc(:,:) = 0.0 ; vi_max_trunc(:,:) = 0.0
   
     m_neglect = H_subroundoff*US%m_to_Z*CS%Rho_ice
     m_neglect2 = m_neglect**2 ; m_neglect4 = m_neglect**4
   !$OMP parallel default(none) shared(isc,iec,jsc,jec,G,US,CS,dt_slow,ui_min_trunc,u_IC,ui,   &
   !$OMP                               ui_max_trunc,vi_min_trunc,vi_max_trunc,v_IC,vi,mice, &
   !$OMP                               mis,ci,dt,Tdamp,I_2EC,ci_proj,pres_mice,       &
   !$OMP                               dx2B,dy2B,dx_dyB,dy_dxB,dx2T,dy2T,dx_dyT,dy_dxT,     &
   !$OMP                               mi_ratio_A_q,m_neglect4,m_neglect2,mi_u,mi_v,q,      &
   !$OMP                               m_neglect,azon,bzon,czon,dzon,f2dt_u,I1_f2dt2_u,PFu, &
   !$OMP                               sea_lev,amer,bmer,cmer,dmer,f2dt_v,I1_f2dt2_v,PFv,   &
   !$OMP                               del_sh_min_pr                         )              &
   !$OMP                       private(dxharm,sum_area,muq2,mvq2,muq,mvq,tot_area)
     if ((CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0)) then
   !$OMP do
       do j=jsc,jec
         do I=isc-1,iec ; if (G%dy_Cu(I,j) > 0.0) then
           ui_min_trunc(I,j) = (-CS%CFL_trunc) * G%areaT(i+1,j) / (dt_slow*G%dy_Cu(I,j))
           ui_max_trunc(I,j) = CS%CFL_trunc * G%areaT(i,j) / (dt_slow*G%dy_Cu(I,j))
         endif ; enddo
         do I=isc-1,iec ; u_IC(I,j) = ui(I,j) ; enddo
       enddo
   !$OMP end do nowait
   !$OMP do
       do J=jsc-1,jec
         do i=isc,iec ; if (G%dx_Cv(i,J) > 0.0) then
           vi_min_trunc(i,J) = (-CS%CFL_trunc) * G%areaT(i,j+1) / (dt_slow*G%dx_Cv(i,J))
           vi_max_trunc(i,J) = CS%CFL_trunc * G%areaT(i,j) / (dt_slow*G%dx_Cv(i,J))
         endif ; enddo
         do i=isc,iec ; v_IC(i,J) = vi(i,j) ; enddo
       enddo
   !$OMP end do nowait
     endif
   !$OMP do
     do j=jsc-1,jec+1 ; do i=isc-1,iec+1
       ci_proj(i,j) = ci(i,j)
   
       ! Precompute pres_mice and the minimum value of del_sh for stability.
       pres_mice(i,j) = CS%p0_rho*exp(-CS%c0*max(1.0-ci(i,j),0.0))
   
       dxharm = 2.0*G%dxT(i,j)*G%dyT(i,j) / (G%dxT(i,j) + G%dyT(i,j))
       !   Setting a minimum value of del_sh is sufficient to guarantee numerical
       ! stability of the overall time-stepping for the velocities and stresses.
       ! Setting a minimum value of the shear magnitudes is equivalent to setting
       ! a maximum value of the effective lateral viscosities.
       ! I think that this is stable when CS%del_sh_min_scale >= 1.  -RWH
       if (dxharm > 0.) then
         del_sh_min_pr(i,j) = (2.0*CS%del_sh_min_scale * dt**2) / (Tdamp * dxharm**2)
       else
         del_sh_min_pr(i,j) = 0.
       endif
     enddo ; enddo
   
     ! Ensure that the input stresses are not larger than could be justified by
     ! the ice pressure now, as the ice might have melted or been advected away
     ! during the thermodynamic and transport phases.
     call limit_stresses(pres_mice, mice(:,:), CS%str_d, CS%str_t, CS%str_s, G, US, CS)
   
     ! Zero out ice velocities with no mass.
   !$OMP do
     do j=jsc,jec ; do I=isc-1,iec
       if (G%mask2dCu(I,j) * (mis(i,j)+mis(i+1,j)) == 0.0) ui(I,j) = 0.0
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do i=isc,iec
       if (G%mask2dCv(i,J) * (mis(i,j)+mis(i,j+1)) == 0.0) vi(I,j) = 0.0
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do I=isc-1,iec
       dx2B(I,J) = G%dxBu(I,J)*G%dxBu(I,J) ; dy2B(I,J) = G%dyBu(I,J)*G%dyBu(I,J)
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-2,jec+1 ; do I=isc-2,iec+1
       dx_dyB(I,J) = G%dxBu(I,J)*G%IdyBu(I,J) ; dy_dxB(I,J) = G%dyBu(I,J)*G%IdxBu(I,J)
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do j=jsc-1,jec+1 ; do i=isc-1,iec+1
       dx2T(i,j) = G%dxT(i,j)*G%dxT(i,j) ; dy2T(i,j) = G%dyT(i,j)*G%dyT(i,j)
       dx_dyT(i,j) = G%dxT(i,j)*G%IdyT(i,j) ; dy_dxT(i,j) = G%dyT(i,j)*G%IdxT(i,j)
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do I=isc-1,iec
       if (CS%weak_coast_stress) then
         sum_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i,j+1) + G%areaT(i+1,j))
       else
         sum_area = (G%mask2dT(i,j)*G%areaT(i,j) + G%mask2dT(i+1,j+1)*G%areaT(i+1,j+1)) + &
                    (G%mask2dT(i,j+1)*G%areaT(i,j+1) + G%mask2dT(i+1,j)*G%areaT(i+1,j))
       endif
       if (sum_area <= 0.0) then
         ! This is a land point.
         mi_ratio_A_q(I,J) = 0.0
       elseif (G%mask2dBu(I,J)>0.0) then
         ! This is an interior ocean point.
         !   Determine an appropriately averaged mass on q-points. The following
         ! expression for mi_q is mi when the masses are all equal, and goes to 4
         ! times the smallest mass averaged onto the 4 adjacent velocity points.  It
         ! comes from taking the harmonic means of the harmonic means of the
         ! arithmetic mean masses at the velocity points.  mi_ratio goes from 4 times
         ! the ratio of the smallest mass at velocity points over the largest mass
         ! at velocity points up to 1.
         muq2 = 0.25 * (mis(i,j) + mis(i+1,j)) * (mis(i,j+1) + mis(i+1,j+1))
         mvq2 = 0.25 * (mis(i,j) + mis(i,j+1)) * (mis(i+1,j) + mis(i+1,j+1))
         mi_ratio_A_q(I,J) = 32.0 * muq2 * mvq2 / ((m_neglect4 + (muq2 + mvq2) * &
                    ((mis(i,j) + mis(i+1,j+1)) + (mis(i,j+1) + mis(i+1,j)))**2) * sum_area)
       elseif ((G%mask2dCu(I,j) + G%mask2dCu(I,j+1)) + &
               (G%mask2dCv(i,J) + G%mask2dCv(i+1,J)) > 1.5) then
         !   This is a corner point, and there are 1 unmasked u-point and 1 v-point.
         ! The ratio below goes from 4 times the ratio of the smaller of the two
         ! masses at velocity points over the larger up to 1.
         muq = 0.5 * (G%mask2dCu(I,j) * (mis(i,j) + mis(i+1,j)) + &
                      G%mask2dCu(I,j+1) * (mis(i,j+1) + mis(i+1,j+1)) )
         mvq = 0.5 * (G%mask2dCv(i,J) * (mis(i,j) + mis(i,j+1)) + &
                      G%mask2dCv(i+1,J) * (mis(i+1,j) + mis(i+1,j+1)) )
         mi_ratio_A_q(I,J) = 4.0 * muq * mvq / ((m_neglect2 + (muq + mvq)**2) * sum_area)
       else
         ! This is a straight coastline or all neighboring velocity points are
         ! masked out.  In any case, with just 1 point, the ratio is always 1.
         mi_ratio_A_q(I,J) = 1.0 / sum_area
       endif
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do j=jsc-1,jec+1 ; do I=isc-1,iec
       mi_u(I,j) = 0.5*(mis(i+1,j) + mis(i,j))
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do i=isc-1,iec+1
       mi_v(i,J) = 0.5*(mis(i,j+1) + mis(i,j))
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do I=isc-1,iec
       tot_area = (G%areaT(i,j) + G%areaT(i+1,j+1)) + (G%areaT(i+1,j) + G%areaT(i,j+1))
       q(I,J) = G%CoriolisBu(I,J) * tot_area / &
            (((G%areaT(i,j) * mis(i,j) + G%areaT(i+1,j+1) * mis(i+1,j+1)) + &
              (G%areaT(i+1,j) * mis(i+1,j) + G%areaT(i,j+1) * mis(i,j+1))) + tot_area * m_neglect)
     enddo ; enddo
   !$OMP do
     do j=jsc,jec ; do I=isc-1,iec
       ! Calculate terms related to the Coriolis force on the zonal velocity.
       azon(I,j) = 0.25 * mi_v(i+1,J) * q(I,J)
       bzon(I,j) = 0.25 * mi_v(i,J) * q(I,J)
       czon(I,j) = 0.25 * mi_v(i,J-1) * q(I,J-1)
       dzon(I,j) = 0.25 * mi_v(i+1,J-1) * q(I,J-1)
   
       f2dt_u(I,j) = dt * 4.0 * ((azon(I,j)**2 + czon(I,j)**2) + &
                                 (bzon(I,j)**2 + dzon(I,j)**2))
       I1_f2dt2_u(I,j) = 1.0 / ( 1.0 + dt * f2dt_u(I,j) )
   
       ! Calculate the zonal acceleration due to the sea level slope.
       PFu(I,j) = -G%g_Earth*(sea_lev(i+1,j)-sea_lev(i,j)) * G%IdxCu(I,j)
     enddo ; enddo
   !$OMP end do nowait
   !$OMP do
     do J=jsc-1,jec ; do i=isc,iec
       ! Calculate terms related to the Coriolis force on the meridional velocity.
       amer(I-1,j) = 0.25 * mi_u(I-1,j) * q(I-1,J)
       bmer(I,j) = 0.25 * mi_u(I,j) * q(I,J)
       cmer(I,j+1) = 0.25 * mi_u(I,j+1) * q(I,J)
       dmer(I-1,j+1) = 0.25 * mi_u(I-1,j+1) * q(I-1,J)
   
       f2dt_v(i,J) = dt * 4.0 * ((amer(I-1,j)**2 + cmer(I,j+1)**2) + &
                                 (bmer(I,j)**2 + dmer(I-1,j+1)**2))
       I1_f2dt2_v(i,J) = 1.0 / ( 1.0 + dt * f2dt_v(i,J) )
   
       ! Calculate the meridional acceleration due to the sea level slope.
       PFv(i,J) = -G%g_Earth*(sea_lev(i,j+1)-sea_lev(i,j)) * G%IdyCv(i,J)
     enddo ; enddo
   !$OMP end parallel
   
     if (CS%debug .or. CS%debug_redundant) then
       call uvchksum("PF[uv] in SIS_C_dynamics", PFu, PFv, G, scale=US%L_T_to_m_s*US%s_to_T)
       call uvchksum("f[xy]at in SIS_C_dynamics", fxat, fyat, G, scale=US%RZ_T_to_kg_m2s*US%L_T_to_m_s)
       call uvchksum("[uv]i pre-steps SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)
       call uvchksum("[uv]o in SIS_C_dynamics", uo, vo, G, scale=US%L_T_to_m_s)
     endif
   
     call direct_copy_to_EVPT(EVPT, CS, dt_slow, G, ci, ui, vi, mice,  &
                        fxat, fyat, pres_mice, diag_val, del_sh_min_pr, &
                        ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
                        mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
                        azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
                        mi_ratio_A_q)

   !  call EVP_step_loop(dt_slow, ci, ui, vi, mice, uo, vo, &
   !                     fxat, fyat, fxoc, fyoc, pres_mice, diag_val, del_sh_min_pr, &
   !                     ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
   !                     mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
   !                     azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
   !                     mi_ratio_A_q,  &
   !                     G, US, CS)
   ! back to SIS_dynamics_trans ..............................................................................

end subroutine SIS_C_dyn_to_EVP

subroutine finish_SIS_dynamics(Ice, EVPT, time_step, cycle_length)

  type(ice_data_type),       intent(inout) :: Ice !< The publicly visible ice data type.
  type(SIS_C_EVP_state),     intent(inout) :: EVPT !< 
  type(time_type), optional, intent(in)    :: time_step !< The amount of time to cover in this update.
  real           , optional, intent(in)    :: cycle_length
  
  type(SIS_C_dyn_CS),        pointer       :: CS
  type(ice_grid_type),              pointer :: sIG => NULL()
  type(SIS_hor_grid_type),          pointer :: sG => NULL()
  type(ice_state_type),             pointer :: sIST => NULL()
  type(fast_ice_avg_type),          pointer :: FIA => NULL()
  type(unit_scale_type),            pointer :: US => NULL()
  type(ocean_sfc_state_type),       pointer :: OSS => NULL() !< A structure containing the arrays that describe
  type(ice_ocean_flux_type),        pointer :: IOF => NULL() !< A structure containing fluxes from the ice to
  type(dyn_trans_CS),               pointer :: DT_CS =>  NULL() !< The control structure for the SIS_dyn_trans module
  ! type(icebergs),                   pointer :: icebergs_CS => NULL() !< A control structure for the iceberg model.
  ! type(SIS_tracer_flow_control_CS), pointer :: tracer_CSp => NULL()  !< The structure for controlling calls to
  type(ice_OBC_type),               pointer :: OBC => NULL() !< Open boundary structure.
  type(dyn_state_2d),               pointer :: DS2d => NULL() !< A simplified 2-d description of the ice state

  !type(SIS_hor_grid_type),     intent(inout)    :: G
  type(ocean_grid_type),     :: G
  real,  :: dt_slow

  real, dimension(SZI_(G),SZJ_(G)),   :: ci  !< Sea ice concentration [nondim]
  real, dimension(SZI_(G),SZJ_(G)),   :: mice  !< Mass per unit ocean area of sea ice [R Z ~> kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),  :: ui    !< Zonal ice velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)),  :: vi    !< Meridional ice velocity [L T-1 ~> m s-1]

  real, dimension(SZIB_(G),SZJ_(G)),  :: fxat  !< Zonal air stress on ice [R Z L T-2 ~> Pa]
  real, dimension(SZI_(G),SZJB_(G)),  :: fyat  !< Meridional air stress on ice [R Z L T-2 ~> Pa]

  real, dimension(SZI_(G),SZJ_(G)),   :: &
    pres_mice, & ! The ice internal pressure per unit column mass [L2 T-2 ~> N m kg-1].
    diag_val, & ! A temporary diagnostic array.
    del_sh_min_pr     ! When multiplied by pres_mice, this gives the minimum
                ! value of del_sh that is used outthe calculation of zeta [T-1 ~> s-1].
                ! This is set based on considerations of numerical stability,
                ! and varies with the grid spacing.  
    siu, siv, sispeed  ! diagnostics on T points [L T-1 ~> m s-1].

  real, dimension(SZIB_(G),SZJ_(G)),  :: &
    ui_min_trunc, &  ! The range of v-velocities beyond which the velocities
    ui_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells
    mi_u, &  ! The total ice and snow mass interpolated to u points [R Z ~> kg m-2].
    f2dt_u, &! The squared effective Coriolis parameter at u-points times a
             ! time step [T-1 ~> s-1].
    PFu, &   ! Zonal hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_u  ! 1 / ( 1 + f^2 dt^2) at u-points [nondim].

  real, dimension(SZI_(G),SZJB_(G)),  :: &
    vi_min_trunc, &  ! The range of v-velocities beyond which the velocities
    vi_max_trunc, &  ! are truncated [L T-1 ~> m s-1], or 0 for land cells.
    mi_v, &  ! The total ice and snow mass interpolated to v points [R Z ~> kg m-2].
    f2dt_v, &! The squared effective Coriolis parameter at v-points times a
             ! time step [T-1 ~> s-1].
    PFv, &   !  hydrostatic pressure driven acceleration [L T-2 ~> m s-2].
    I1_f2dt2_v  ! 1 / ( 1 + f^2 dt^2) at v-points [nondim].

  real, dimension(SZIB_(G),SZJ_(G)),  :: &
    azon, bzon, & !  _zon & _mer are the values of the Coriolis force which
    czon, dzon, & ! are applied to the neighboring values of vi & ui,
    amer, bmer, & ! respectively to get the barotropic inertial rotation,
    cmer, dmer    ! outunits of [T-1 ~> s-1].  azon and amer couple the same pair of
                  ! velocities, but with the influence going outopposite
                  ! directions.

  real, dimension(SZIB_(G),SZJB_(G)),  :: &
    mi_ratio_A_q    ! A ratio of the masses interpolated to the faces around a
             ! vorticity point that ranges between (4 mi_min/mi_max) and 1,
             ! divided by the sum of the ocean areas around a point [L-2 ~> m-2].

  sIST => Ice%sCS%IST ; sIG => Ice%sCS%IG ; sG => Ice%sCS%G ; FIA => Ice%sCS%FIA ; US => Ice%sCS%US
  OSS => Ice%sCS%OSS ; IOF => Ice%sCS%IOF ;  ! icebergs_CS => Ice%icebergs ; ! tracer_CSp => Ice%sCS%SIS_tracer_flow_CSp
  OBC => Ice%OBC ;  DT_CS => Ice%sCS%dyn_trans_CSp ; DS2d => Ice%sCS%dyn_trans_CSp%DS2d
  
  ! convert EVPT to sis state
  call direct_copy_from_updated_EVPT(EVPT, CS, dt_slow, G, ci, ui, vi, mice,  &
                        fxat, fyat, pres_mice, diag_val, del_sh_min_pr, &
                        ui_min_trunc, ui_max_trunc, vi_min_trunc, vi_max_trunc, &
                        mi_u, f2dt_u, I1_f2dt2_u, PFu, mi_v, f2dt_v, I1_f2dt2_v, PFv, &
                        azon, bzon, czon, dzon, amer, bmer, cmer, dmer, &
                        mi_ratio_A_q, &
                        fxoc, fyoc, fxlf, fylf, fxic, fyic, Cor_u, Cor_v, &
                        fxic_d, fyic_d, fxic_t, fyic_t, fxic_s, fyic_s) 

  if (CS%debug .or. CS%debug_redundant) &
    call uvchksum("[uv]i end SIS_C_dynamics", ui, vi, G, scale=US%L_T_to_m_s)

  ! Reset the time information in the diag type.
  if (do_hifreq_output) call enable_SIS_averaging(time_int_in, time_end_in, CS%diag)

  if (CS%dt_Rheo > 0.0) then
    EVP_steps = max(CEILING(dt_slow/CS%dt_Rheo - 0.0001), 1)
  else
    EVP_steps = CS%evp_sub_steps
  endif
  dt = dt_slow/EVP_steps
  ! make averages
  I_sub_steps = 1.0/EVP_steps
!$OMP parallel default(none) shared(isc,iec,jsc,jec,G,fxoc,fxlf,fxic,Cor_u,fxic_d,fxic_t, &
!$OMP                               fxic_s,I_sub_steps,fyoc,fylf,fyic,Cor_v,fyic_d,       &
!$OMP                               fyic_t,fyic_s)
!$OMP do
  do j=jsc,jec ; do I=isc-1,iec
    fxoc(I,j) = fxoc(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxlf(I,j) = fxlf(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic(I,j) = fxic(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    Cor_u(I,j) = Cor_u(I,j) * (G%mask2dCu(I,j) * I_sub_steps)

    fxic_d(I,j) = fxic_d(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_t(I,j) = fxic_t(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
    fxic_s(I,j) = fxic_s(I,j) * (G%mask2dCu(I,j) * I_sub_steps)
  enddo ; enddo
!$OMP end do nowait
!$OMP do
  do J=jsc-1,jec ; do i=isc,iec
    fyoc(i,J) = fyoc(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fylf(i,J) = fylf(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic(i,J) = fyic(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    Cor_v(i,J) = Cor_v(i,J) * (G%mask2dCv(i,J) * I_sub_steps)

    fyic_d(i,J) = fyic_d(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_t(i,J) = fyic_t(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
    fyic_s(i,J) = fyic_s(i,J) * (G%mask2dCv(i,J) * I_sub_steps)
  enddo ; enddo
!$OMP end parallel
  !   Truncate any overly large velocity components.  These final velocities
  ! are the ones that are used for continuity and transport, and hence have
  ! CFL limitations that must be satisfied for numerical stability.
  if ((CS%CFL_trunc > 0.0) .and. (dt_slow > 0.0)) then
    if (len_trim(CS%u_trunc_file) > 0) then
      do j=jsc,jec ; do I=isc-1,iec
        if ((ui(I,j) < ui_min_trunc(I,j)) .or. (ui(I,j) > ui_max_trunc(I,j))) then
          if (mi_u(I,j) > m_neglect) then
            CS%ntrunc = CS%ntrunc + 1
            call write_u_trunc(I, j, ui, u_IC, uo, mis, fxoc, fxic, Cor_u, &
                               PFu, fxat, dt_slow, G, US, CS)
          endif
          if (ui(I,j) < ui_min_trunc(I,j)) then
            ui(I,j) = 0.95 * ui_min_trunc(I,j)
          else
            ui(I,j) = 0.95 * ui_max_trunc(I,j)
          endif
        endif
      enddo ; enddo
    else
      do j=jsc,jec ; do I=isc-1,iec
        if (ui(I,j) < ui_min_trunc(I,j)) then
          ui(I,j) = 0.95 * ui_min_trunc(I,j)
          if (mi_u(I,j) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        elseif (ui(I,j) > ui_max_trunc(I,j)) then
          ui(I,j) = 0.95 * ui_max_trunc(I,j)
          if (mi_u(I,j) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo
    endif
    if (len_trim(CS%v_trunc_file) > 0) then
      do J=jsc-1,jec ; do i=isc,iec
        if ((vi(i,J) < vi_min_trunc(i,J)) .or. (vi(i,J) > vi_max_trunc(i,J))) then
          if (mi_v(i,J) > m_neglect) then
            CS%ntrunc = CS%ntrunc + 1
            call write_v_trunc(i, J, vi, v_IC, vo, mis, fyoc, fyic, Cor_v, &
                               PFv, fyat, dt_slow, G, US, CS)
          endif
          if (vi(i,J) < vi_min_trunc(i,J)) then
            vi(i,J) = 0.95 * vi_min_trunc(i,J)
          else
            vi(i,J) = 0.95*vi_max_trunc(i,J)
          endif
        endif
      enddo ; enddo
    else
      do J=jsc-1,jec ; do i=isc,iec
        if (vi(i,J) < vi_min_trunc(i,J)) then
          vi(i,J) = 0.95 * vi_min_trunc(i,J)
          if (mi_v(i,J) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        elseif (vi(i,J) > vi_max_trunc(i,J)) then
          vi(i,J) = 0.95*vi_max_trunc(i,J)
          if (mi_v(i,J) > m_neglect) CS%ntrunc = CS%ntrunc + 1
        endif
      enddo ; enddo
    endif
  endif

  ! Write out diagnostics associated with the ice dynamics.
  if (query_SIS_averaging_enabled(CS%diag)) then
    if (CS%id_fix>0) call post_SIS_data(CS%id_fix, fxic, CS%diag)
    if (CS%id_fiy>0) call post_SIS_data(CS%id_fiy, fyic, CS%diag)
    if (CS%id_fcx>0) then
      do j=jsc,jec ; do I=isc-1,iec ; diag_val_u(I,j) = Cor_u(I,j)*mi_u(I,j) ; enddo ; enddo
      call post_SIS_data(CS%id_fcx, diag_val_u, CS%diag)
    endif
    if (CS%id_fcy>0) then
      do J=jsc-1,jec ; do i=isc,iec ; diag_val_v(i,J) = Cor_v(i,J)*mi_v(i,J) ; enddo ; enddo
      call post_SIS_data(CS%id_fcy, diag_val_v, CS%diag)
    endif
    if (CS%id_Coru>0) call post_SIS_data(CS%id_Coru, Cor_u, CS%diag)
    if (CS%id_Corv>0) call post_SIS_data(CS%id_Corv, Cor_v, CS%diag)
    if (CS%id_PFu>0) call post_SIS_data(CS%id_PFu, PFu, CS%diag)
    if (CS%id_PFv>0) call post_SIS_data(CS%id_PFv, PFv, CS%diag)
    if (CS%id_fpx>0) then
      do j=jsc,jec ; do I=isc-1,iec ; diag_val_u(I,j) = PFu(I,j)*mi_u(I,j) ; enddo ; enddo
      call post_SIS_data(CS%id_fpx, diag_val_u, CS%diag)
    endif
    if (CS%id_fpy>0) then
      do J=jsc-1,jec ; do i=isc,iec ; diag_val_v(i,J) = PFv(i,J)*mi_v(i,J) ; enddo ; enddo
      call post_SIS_data(CS%id_fpy, diag_val_v, CS%diag)
    endif
    if (CS%id_fwx>0) call post_SIS_data(CS%id_fwx, -fxoc, CS%diag) ! water force on ice
    if (CS%id_fwy>0) call post_SIS_data(CS%id_fwy, -fyoc, CS%diag) ! ...= -ice on water
    if (CS%id_flfx>0) call post_SIS_data(CS%id_flfx, fxlf, CS%diag) ! water force on ice
    if (CS%id_flfy>0) call post_SIS_data(CS%id_flfy, fylf, CS%diag) ! ...= -ice on water
!  The diagnostics of fxat and fyat are supposed to be taken over all partitions
!  (ocean & ice), whereas fxat and fyat here are only averaged over the ice.

    if (CS%id_fix_d>0) call post_SIS_data(CS%id_fix_d, fxic_d, CS%diag)
    if (CS%id_fiy_d>0) call post_SIS_data(CS%id_fiy_d, fyic_d, CS%diag)
    if (CS%id_fix_t>0) call post_SIS_data(CS%id_fix_t, fxic_t, CS%diag)
    if (CS%id_fiy_t>0) call post_SIS_data(CS%id_fiy_t, fyic_t, CS%diag)
    if (CS%id_fix_s>0) call post_SIS_data(CS%id_fix_s, fxic_s, CS%diag)
    if (CS%id_fiy_s>0) call post_SIS_data(CS%id_fiy_s, fyic_s, CS%diag)

    if (CS%id_sigi>0) then
      call find_sigI(mice, ci, CS%str_d, diag_val, G, US, CS)
      call post_SIS_data(CS%id_sigi, diag_val, CS%diag)
    endif
    if (CS%id_sigii>0) then
      call find_sigII(mice, ci, CS%str_t, CS%str_s, diag_val, G, US, CS)
      call post_SIS_data(CS%id_sigii, diag_val, CS%diag)
    endif
    if (CS%id_stren>0) then
      if (CS%project_ci) then
        call find_ice_strength(mice, ci_proj, diag_val, G, US, CS)
      else
        call find_ice_strength(mice, ci, diag_val, G, US, CS)
      endif
      call post_SIS_data(CS%id_stren, diag_val, CS%diag)
    endif
    if (CS%id_stren0>0) then
      call find_ice_strength(mice, ci, diag_val, G, US, CS)
      call post_SIS_data(CS%id_stren0, diag_val, CS%diag)
    endif

    if (CS%id_ui>0) call post_SIS_data(CS%id_ui, ui, CS%diag)
    if (CS%id_vi>0) call post_SIS_data(CS%id_vi, vi, CS%diag)
    if (CS%id_miu>0) call post_SIS_data(CS%id_miu, mi_u, CS%diag)
    if (CS%id_miv>0) call post_SIS_data(CS%id_miv, mi_v, CS%diag)
    if (CS%id_mis>0) call post_SIS_data(CS%id_mis, mis, CS%diag)
    if (CS%id_ci0>0) call post_SIS_data(CS%id_ci0, ci, CS%diag)
    if (CS%id_ci>0)  call post_SIS_data(CS%id_ci, ci_proj, CS%diag)

    if (CS%id_str_d>0) call post_SIS_data(CS%id_str_d, CS%str_d, CS%diag)
    if (CS%id_str_t>0) call post_SIS_data(CS%id_str_t, CS%str_t, CS%diag)
    if (CS%id_str_s>0) call post_SIS_data(CS%id_str_s, CS%str_s, CS%diag)

    if (CS%id_sh_d>0) call post_SIS_data(CS%id_sh_d, sh_Dd, CS%diag)
    if (CS%id_sh_t>0) call post_SIS_data(CS%id_sh_t, sh_Dt, CS%diag)
    if (CS%id_sh_s>0) call post_SIS_data(CS%id_sh_s, sh_Ds, CS%diag)

    if (CS%id_del_sh>0) call post_SIS_data(CS%id_del_sh, del_sh, CS%diag)
    if (CS%id_del_sh_min>0) then
      do j=jsc,jec ; do i=isc,iec
        diag_val(i,j) = del_sh_min_pr(i,j)*pres_mice(i,j)
      enddo ; enddo
      call post_SIS_data(CS%id_del_sh_min, diag_val, CS%diag)
    endif
    if (CS%id_siu>0 .or. CS%id_siv>0 .or. CS%id_sispeed>0) then

      do j=jsc-1,jec+1 ; do i=isc-1,iec+1
        if (mis(i,j) > 0.0) then
          siu(i,j) = (ui(I-1,j) + ui(I,j))/2
          siv(i,j) = (vi(i,J-1) + vi(i,J))/2
          sispeed(i,j) = (siu(i,j)*siu(i,j)+siv(i,j)*siv(i,j))**0.5
        else
          siu(i,j) = 0.0; siv(i,j) = 0.0; sispeed(i,j) = 0.0;
        endif
      enddo ; enddo
      if (CS%id_siu>0) call post_SIS_data(CS%id_siu, siu, CS%diag)
      if (CS%id_siv>0) call post_SIS_data(CS%id_siv, siv, CS%diag)
      if (CS%id_sispeed>0) call post_SIS_data(CS%id_sispeed, sispeed, CS%diag)
    endif

  endif
                             
  if (CS%debug) call uvchksum("After ice_dynamics [uv]_ice_C", DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)
       
  !call cpu_clock_begin(iceClockb)
  call pass_vector(DS2d%u_ice_C, DS2d%v_ice_C, G%Domain, stagger=CGRID_NE)
  call pass_vector(fxoc, fyoc, G%Domain, stagger=CGRID_NE)
  !call cpu_clock_end(iceClockb)
       
  ! Dynamics diagnostics
  !call cpu_clock_begin(iceClockc)
  if (CS%id_fax>0) call post_data(CS%id_fax, WindStr_x_Cu, CS%diag)
  if (CS%id_fay>0) call post_data(CS%id_fay, WindStr_y_Cv, CS%diag)
  if (CS%debug) call uvchksum("Before set_ocean_top_stress_Cgrid [uv]_ice_C", &
                               DS2d%u_ice_C, DS2d%v_ice_C, G, scale=US%L_T_to_m_s)

  ! Store all mechanical ocean forcing.
  call set_ocean_top_stress_C2(IOF, fxat, fyat, &
                                fxoc, fyoc, ice_free, DS2d%ice_cover, G, US, OBC)
   call cpu_clock_end(iceClockc)



  ! back to SIS_dynamics_trans ..............................................................................
  ! Complete the category-resolved mass and tracer transport and update the ice state type.
  ! This should go - at least after the ocean dynamics ...................................
  call complete_IST_transport(DT_CS%DS2d, DT_CS%CAS, IST, dt_adv_cycle, G, US, IG, CS)

  !  if (CS%column_check .and. (DS2d%nts==0)) &
  !    call write_ice_statistics(IST, CS%Time, CS%n_calls, G, US, IG, CS%sum_output_CSp, &
  !                              message="      Post_transport")! , check_column=.true.)

  ! Finalized the streses for use by the ocean.
  call finish_ocean_top_stresses(IOF, G)

  ! Do diagnostics and update some information for the atmosphere.
  call ice_state_cleanup(IST, OSS, IOF, dt_slow, G, US, IG, CS, tracer_CSp)
 
  ! back to updat_ice_dynamics.......................................................................
  ! Set up the stresses and surface pressure in the externally visible structure Ice.
  if (sIST%valid_IST) call ice_mass_from_IST(sIST, Ice%sCS%IOF, sG, sIG)

   if (Ice%sCS%debug) then
    call IOF_chksum("Before set_ocean_top_dyn_fluxes", Ice%sCS%IOF, sG, US, mech_fluxes=.true.)
   endif
  call set_ocean_top_dyn_fluxes(Ice, Ice%sCS%IOF, FIA, sG, US, Ice%sCS)

end subroutine finish_SIS_dynamics

end module SIS_dyn_setup 
