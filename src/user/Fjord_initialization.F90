module Fjord_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_ALE_sponge, only : ALE_sponge_CS, set_up_ALE_sponge_field, initialize_ALE_sponge
use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_io, only : file_exists
use MOM_io, only : MOM_read_data
use MOM_io, only : slasher
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA
use regrid_consts, only : REGRIDDING_SIGMA_SHELF_ZSTAR

use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_SIMPLE
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_tracer_registry, only : tracer_name_lookup



implicit none ; private

#include <MOM_memory.h>

! -----------------------------------------------------------------------------
! Private (module-wise) parameters
! -----------------------------------------------------------------------------

character(len=40) :: mdl = "Fjord_initialization" ! This module's name.

! -----------------------------------------------------------------------------
! The following routines are visible to the outside world
! -----------------------------------------------------------------------------
public Fjord_initialize_topography
public Fjord_initialize_thickness
public Fjord_initialize_temperature_salinity
public Fjord_initialize_sponges
public Fjord_set_OBC_data
! -----------------------------------------------------------------------------
! This module contains the following routines
! -----------------------------------------------------------------------------
contains

!> Initialization of topography
subroutine Fjord_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),             intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                      intent(out) :: D !< Ocean bottom depth in m
  type(param_file_type),              intent(in)  :: param_file !< Parameter file structure
  real,                               intent(in)  :: max_depth  !< Maximum depth of model in m

! This subroutine sets up the ISOMIP topography
  real :: min_depth ! The minimum and maximum depths in m.

! The following variables are used to set up the bathymetry in the ISOMIP example.
! check this paper: http://www.geosci-model-dev-discuss.net/8/9859/2015/gmdd-8-9859-2015.pdf
!GM ISOMIP:
! real :: bmax            ! max depth of bedrock topography
! real :: b0,b2,b4,b6     ! first, second, third and fourth bedrock topography coeff
! real :: xbar           ! characteristic along-flow lenght scale of the bedrock
! real :: dc              ! depth of the trough compared with side walls
! real :: fc              ! characteristic width of the side walls of the channel
! real :: wc              ! half-width of the trough
! real :: ly              ! domain width (across ice flow)
! real :: bx, by, xtil    ! dummy vatiables
 logical :: is_2D         ! If true, use 2D setup


! TJM Fjord topo varibles

  real :: fw            ! width of the fjord in the i,x,zonal coordinate
  real :: fl            ! length of fjord in the j,y,meridional coordinate
  real :: shl           ! length of shelf in j
  real :: sll           ! length of the shelf slope in j
  real :: osl           ! length of off shelf region
  real :: shw           ! width of shelf reginon in i (and entire domain) 

  real :: osd           ! depth at the botom of the slope
  real :: shd           ! depth on the shelf
  real :: fd            ! depth in the fjord

  real :: land=0.0      ! where there is land
  real :: fs, fe        ! start/end position of the fjord in x        

  real :: m, b          ! slop parametres

  real :: xtil,ytil

! G%ieg and G%jeg are the last indices in the global domain

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "Fjord_initialize_topography" ! This subroutine's name.
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed


  call MOM_mesg(" Fjord_initialization.F90, Fjord_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)
  call get_param(param_file, mdl, "Fjord_2D",is_2D,'If true, use a 2D setup.', default=.false.)

!TJM Read topo vars from param_file
   
  call get_param(param_file, mdl, "Fjord_width", fw, &
                 "The width of the fjord", units="m", default=10.0E3)
  call get_param(param_file, mdl, "Fjord_length", fl, &
                 "The length of the fjord", units="m", default=60.0E3)
  call get_param(param_file, mdl, "Shelf_length", shl, &
                 "The length of the Shelf", units="m", default=10.0E3)
  call get_param(param_file, mdl, "Slope_length", sll, &
                 "The length of the continental shelf slope", units="m", default=10.0E3)
  call get_param(param_file, mdl, "OffShelf_length", osl, &
                 "The length of the off shelf region", units="m", default=10.0E3)
  call get_param(param_file, mdl, "Shelf_width", shw, &
                 "The width of the shelf", units="m", default=40.0E3)

  call get_param(param_file, mdl, "Off_Shelf_depth", osd, &
                 "The depth off the shelf", units="m", default=2000.0)
  call get_param(param_file, mdl, "Shelf_depth", shd, &
                 "The depth of the shelf", units="m", default=500.0)
  call get_param(param_file, mdl, "Fjord_depth", fd, &
                 "The depth of the fjord", units="m", default=500.0)

  call get_param(param_file, mdl, "ISOMIP_2D",is_2D,'If true, use a 2D setup.', default=.false.)


!---TJM

!gm The following variables should be transformed into runtime parameters?
!  bmax=720.0; b0=-150.0; b2=-728.8; b4=343.91; b6=-50.57
!  xbar=300.0E3; dc=500.0; fc=4.0E3; wc=24.0E3; ly=80.0E3
!  bx = 0.0; by = 0.0; xtil = 0.0


fs = (shw-fw)/2.
fe = fs + fw

m = (shd - osd)/sll
b = osd

land =-1.0

  if (is_2D) then
    do j=js,je ; do i=is,ie
      ! 2D setup
        ytil = G%geoLatT(i,j)*1.0e3
      if  (ytil<fl) then
        D(i,j)=fd
      elseif (ytil>=fl .and. ytil<shl+fl) then
        D(i,j)=shd
      elseif (ytil>=shl+fl .and. ytil<shl+fl+sll) then
        D(i,j)=m*ytil+b
      elseif (ytil>=fl+shl+sll) then
        D(i,j)=osd
      endif
    enddo ; enddo

  else
    do j=js,je ; do i=is,ie
        xtil = G%geoLonT(i,j)*1.0e3
        ytil = G%geoLatT(i,j)*1.0e3
      if  (ytil>osl+sll+shl) then
        if (xtil<=fs) D(i,j)=land
        if (xtil>fs .and. xtil<fe) D(i,j)=fd
        if (xtil>=fe) D(i,j)=land
!         D(i,j)=fd
      elseif (ytil>=osl+sll .and. ytil<osl+sll+shl) then
        D(i,j)=shd
      elseif (ytil>=osl .and. ytil<osl+sll) then
        D(i,j)=m*(ytil-osl)+b
      elseif (ytil<osl) then
        D(i,j)=osd
      endif

!   if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
    if (D(i,j) > max_depth) D(i,j) = max_depth
    if (D(i,j) < min_depth) D(i,j) = 0.5*min_depth
    
    enddo ; enddo
  endif

end subroutine Fjord_initialize_topography
! -----------------------------------------------------------------------------

!> Initialization of thicknesses
subroutine Fjord_initialize_thickness ( h, G, GV, param_file, tv, just_read_params)
  type(ocean_grid_type),   intent(in)  :: G           !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV          !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized, in H.
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  type(thermo_var_ptrs),   intent(in)  :: tv          !< A structure containing pointers to any
                                                      !! available thermodynamic fields, including
                                                      !! the eqn. of state.
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing h.

  ! Local variables
  real :: e0(SZK_(G)+1)     ! The resting interface heights, in m, usually !
                          ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)! Interface height relative to the sea surface !
                          ! positive upward, in m.                       !
  integer :: i, j, k, is, ie, js, je, nz, tmp1
  real    :: x
  real    :: delta_h, rho_range
  real    :: min_thickness, s_sur, s_bot, t_sur, t_bot, rho_sur, rho_bot
  logical :: just_read    ! If true, just read parameters but set nothing.
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  if (.not.just_read) &
    call MOM_mesg("MOM_initialization.F90, initialize_thickness_uniform: setting thickness")

  call get_param(param_file, mdl,"MIN_THICKNESS",min_thickness, &
                 'Minimum layer thickness', units='m', default=1.e-3, do_not_log=just_read)
  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)

  select case ( coordinateMode(verticalCoordinate) )

!  case ( REGRIDDING_LAYER, REGRIDDING_RHO ) ! Initial thicknesses for isopycnal coordinates

  case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR )   ! Initial thicknesses for z coordinates
    if (just_read) return ! All run-time parameters have been read, so return.
    do j=js,je ; do i=is,ie
      eta1D(nz+1) = -1.0*G%bathyT(i,j)
      do k=nz,1,-1
        eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
        if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
          eta1D(k) = eta1D(k+1) + min_thickness
          h(i,j,k) = GV%m_to_H * min_thickness
        else
          h(i,j,k) = GV%m_to_H * (eta1D(k) - eta1D(k+1))
        endif
      enddo
   enddo ; enddo

 ! case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates
  case default    
    call MOM_error(FATAL,"Fjord_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

end subroutine Fjord_initialize_thickness

!> Initial values for temperature and salinity
subroutine Fjord_initialize_temperature_salinity ( T, S, h, G, GV, param_file, &
                                                    eqn_of_state, just_read_params)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< Vertical grid structure
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: T !< Potential temperature (degC)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(out) :: S !< Salinity (ppt)
  real, dimension(SZI_(G),SZJ_(G), SZK_(G)), intent(in)  :: h !< Layer thickness in H (m or kg m-2)
  type(param_file_type),                     intent(in)  :: param_file !< Parameter file structure
  type(EOS_type),                            pointer     :: eqn_of_state !< Equation of state structure
  logical,       optional, intent(in)  :: just_read_params !< If present and true, this call will
                                                      !! only read parameters without changing T & S.
  character(len=40)  :: mdl = "Fjord_initialize_temperatue_salinity" !This subroutine's name.

  ! Local variables
  integer   :: i, j, k, is, ie, js, je, nz, itt
  real      :: z,z0          ! vertical position in z space
  logical :: just_read    ! If true, just read parameters but set nothing.
  real :: sl_sbu, sl_smu, sl_sbl, sl_sml, sl_tbu, sl_tmu, sl_tbl, sl_tml ! slopes and linear stuff
  real, dimension(SZK_(G)) :: Sf, Tf, Ssl, Tsl, Sos, Tos ! place holders 
  character(len=40) :: verticalCoordinate
  real :: xtil,ytil

  ! Read in vars
  ! Other
  real :: shl           ! length of shelf in j
  real :: sll           ! length of the shelf slope in j
  real :: osl           ! length of off shelf region
  real :: osd           ! depth at the botom of the slope
  real :: shw, fw, fs, fe, fd
  ! Fjord
  real :: ft_sur, fs_sur
  real :: fs_L, ft_L1, ft_L2, ft_sc1, ft_sc2, fs_sc
  real :: fs_m, fj_zc
  real :: s1, s2, t1, t2, buff, zc

  ! Slope
  real :: sl_zc1, sl_zc2
  real :: sl_S0, sl_SD, sl_Smax, sl_T0, sl_tcl, sl_TD
  real :: zc1, zc2 

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  just_read = .false. ; if (present(just_read_params)) just_read = just_read_params

  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
!  call get_param(param_file, mdl, "Fjord_T_SUR",t_sur, &
!                 'Temperature at the surface (interface)', default=-1.9, do_not_log=just_read)

! Other
call get_param(param_file, mdl, "Off_Shelf_depth", osd, &
                 "The depth off the shelf", units="m", default=2000.0)
  call get_param(param_file, mdl, "Slope_length", sll, &
                 "The length of the continental shelf slope", units="m", default=10.0E3)
  call get_param(param_file, mdl, "OffShelf_length", osl, &
                 "The length of the off shelf region", units="m", default=10.0E3)
call get_param(param_file, mdl, "Fjord_width", fw, &
                 "The width of the fjord", units="m", default=10.0E3)
call get_param(param_file, mdl, "Shelf_width", shw, &
                 "The width of the shelf", units="m", default=30.0E3) 
call get_param(param_file, mdl, "Fjord_depth", fd, &
                 "The depth of the fjord", units="m", default=500.0)

call get_param(param_file, mdl, "fj_t1", t1, &
                 "fjord surf temp", default=-1.0, do_not_log=just_read)
call get_param(param_file, mdl, "fj_t2", t2, &
                 "fjord bot temp", default=4.0, do_not_log=just_read)
call get_param(param_file, mdl, "fj_s1", s1, &
                 "fjord surf salt", default=32.0, do_not_log=just_read)
call get_param(param_file, mdl, "fj_s2", s2, &
                 "fjord bot salt", default=35.0, do_not_log=just_read)
call get_param(param_file, mdl, "fj_buff", buff, &
                 "fjord buffer height", default=100.0, do_not_log=just_read)
call get_param(param_file, mdl, "fj_zc", zc, &
                 "fjord interface depth", default=125.0, do_not_log=just_read)

  select case ( coordinateMode(verticalCoordinate) )

   case (  REGRIDDING_RHO, REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR, REGRIDDING_SIGMA )

   fs = shw/2.-fw/2.
   fe = fs + fw

   do j=js,je ; do i=is,ie
        xtil = G%geoLonT(i,j)*1.0e3
        ytil = G%geoLatT(i,j)*1.0e3
        z0=0.0
        do k=1,nz
!          z = z0 + 0.5*h(i,j,k) * GV%H_to_m
!          if (z<zc-buff) then
!              T(i,j,k) = t1
!              S(i,j,k) = s1
!          elseif (z>=zc-buff .and. z<=zc+buff) then
!              T(i,j,k) = (t2-t1)/(2*buff)*(z-(zc-buff))+t1
!              S(i,j,k) = (s2-s1)/(2*buff)*(z-(zc-buff))+s1
!          elseif (z>zc+buff) then
!              T(i,j,k) = t2
!              S(i,j,k) = s2
!          endif
!          z0 = z0 +  h(i,j,k) * GV%H_to_m
!      enddo

        z = z0 + 0.5*h(i,j,k) * GV%H_to_m
      
          if (z>=250) then
             T(i,j,k) = 3.46      
            S(i,j,k) = 34.7
          else
             T(i,j,k) = 2.43
             S(i,j,k) = 34.35
          endif    

         z0 = z0 +  h(i,j,k) * GV%H_to_m
       enddo
      
     enddo ; enddo


!    S(:,:,:) = 35;
!    T(:,:,:) = 1;


!    case ( REGRIDDING_LAYER )

   case default
      call MOM_error(FATAL,"isomip_initialize: "// &
      "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

  end select

  ! for debugging
  !i=G%iec; j=G%jec
  !do k = 1,nz
  !   call calculate_density(T(i,j,k),S(i,j,k),0.0,rho_tmp,eqn_of_state)
  !   write(*,*) 'k,h,T,S,rho,Rlay',k,h(i,j,k),T(i,j,k),S(i,j,k),rho_tmp,GV%Rlay(k)
  !enddo

end subroutine Fjord_initialize_temperature_salinity


!> Sets up the the inverse restoration time (Idamp), and
! the values towards which the interface heights and an arbitrary
! number of tracers should be restored within each sponge.
subroutine Fjord_initialize_sponges(G, GV, tv, PF, use_ALE, CSp, ACSp)
  type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV !< The ocean's vertical grid structure.
  type(thermo_var_ptrs), intent(in) :: tv   !< A structure containing pointers
                                            !! to any available thermodynamic
                                            !! fields, potential temperature and
                                            !! salinity or mixed layer density.
                                            !! Absent fields have NULL ptrs.
  type(param_file_type), intent(in) :: PF   !< A structure indicating the
                                            !! open file to parse for model
                                            !! parameter values.
  logical, intent(in) :: use_ALE            !< If true, indicates model is in ALE mode
  type(sponge_CS),   pointer    :: CSp      !< Layer-mode sponge structure
  type(ALE_sponge_CS),   pointer    :: ACSp !< ALE-mode sponge structure

! Local variables
  real :: T(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for temp
  real :: S(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for salt
  real :: RHO(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for RHO
  real :: h(SZI_(G),SZJ_(G),SZK_(G))  ! A temporary array for thickness
  real :: Idamp(SZI_(G),SZJ_(G))    ! The inverse damping rate, in s-1.
  real :: TNUDG                     ! Nudging time scale, days
  real      :: S_sur, T_sur;        ! Surface salinity and temerature in sponge
  real      :: S_bot, T_bot;        ! Bottom salinity and temerature in sponge
  real      :: t_ref, s_ref         ! reference T and S
  real      :: rho_sur, rho_bot, rho_range, t_range, s_range

  real :: e0(SZK_(G)+1)               ! The resting interface heights, in m, usually !
                                    ! negative because it is positive upward.      !
  real :: eta1D(SZK_(G)+1)          ! Interface height relative to the sea surface !
  real :: eta(SZI_(G),SZJ_(G),SZK_(G)+1) ! A temporary array for eta.

                                    ! positive upward, in m.
  real :: min_depth, dummy1, z, delta_h
  real :: damp, rho_dummy, min_thickness, rho_tmp, xi0
  character(len=40) :: verticalCoordinate, filename, state_file
  character(len=40) :: temp_var, salt_var, eta_var, inputdir

  character(len=40)  :: mdl = "Fjord_initialize_sponges" ! This subroutine's name.

  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz

!--TJM
  logical :: just_read
  real :: sp_wid
  real :: z0 
  real :: fs_sur, fs_L, fs_sc, fs_m, fd
  real :: s1, s2, t1, t2, buff, zc
  real :: xtil, ytil
  real :: shw, shl, sll, osl, land
  logical :: REENTRANT_X

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call get_param(PF, mdl,"MIN_THICKNESS",min_thickness,'Minimum layer thickness',units='m',default=1.e-3)

  call get_param(PF, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
            default=DEFAULT_COORDINATE_MODE)

  call get_param(PF, mdl, "Fjord_TNUDG", TNUDG, 'Nudging time scale for sponge layers (days)',  default=0.0)

  call get_param(PF, mdl, "T_REF", t_ref, 'Reference temperature',  default=10.0,&
                 do_not_log=.true.)

  call get_param(PF, mdl, "S_REF", s_ref, 'Reference salinity',  default=35.0,&
                 do_not_log=.true.)

!--TJM
  call get_param(PF, mdl, "Fjord_WID_SPONGE", sp_wid, 'width of sponge layer.',  default=100.0)


call get_param(PF, mdl, "fj_t1", t1, &
                 "fjord surf temp", default=-1.0, do_not_log=just_read)
call get_param(PF, mdl, "fj_t2", t2, &
                 "fjord bot temp", default=4.0, do_not_log=just_read)
call get_param(PF, mdl, "fj_s1", s1, &
                 "fjord surf salt", default=32.0, do_not_log=just_read)
call get_param(PF, mdl, "fj_s2", s2, &
                 "fjord bot salt", default=35.0, do_not_log=just_read)
call get_param(PF, mdl, "fj_buff", buff, &
                 "fjord buffer height", default=100.0, do_not_log=just_read)
call get_param(PF, mdl, "fj_zc", zc, &
                 "fjord interface depth", default=125.0, do_not_log=just_read)

call get_param(PF, mdl, "Shelf_length", shl, &
                 "The length of the Shelf", units="m", default=10.0E3)
  call get_param(PF, mdl, "Slope_length", sll, &
                 "The length of the continental shelf slope", units="m", default=10.0E3)
  call get_param(PF, mdl, "OffShelf_length", osl, &
                 "The length of the off shelf region", units="m", default=10.0E3)
  call get_param(PF, mdl, "Shelf_width", shw, &
                 "The width of the shelf", units="m", default=30.0E3)

call get_param(PF, mdl, "REENTRANT_X", REENTRANT_X)

  T(:,:,:) = 0.0 ; S(:,:,:) = 0.0 ; Idamp(:,:) = 0.0; RHO(:,:,:) = 0.0
  S_range = s_sur - s_bot
  T_range = t_sur - t_bot

!   Set up sponges for ISOMIP configuration
  call get_param(PF, mdl, "MINIMUM_DEPTH", min_depth, &
                 "The minimum depth of the ocean.", units="m", default=0.0)

   if (associated(CSp)) call MOM_error(FATAL, &
          "Fjord_initialize_sponges called with an associated control structure.")
   if (associated(ACSp)) call MOM_error(FATAL, &
          "Fjord_initialize_sponges called with an associated ALE-sponge control structure.")

  !  Here the inverse damping time, in s-1, is set. Set Idamp to 0     !
  !  wherever there is no sponge, and the subroutines that are called  !
  !  will automatically set up the sponges only where Idamp is positive!
  !  and mask2dT is 1.

   land = osl + sll + shl 

   do i=is,ie; do j=js,je
      xtil = G%geoLonT(i,j)*1.0e3
      ytil = G%geoLatT(i,j)*1.0e3
!      if (REENTRANT_X) then
!            if (0<=ytil .AND. ytil<=sp_wid) then
!
!            ! 1 / day
!            dummy1=(xtil)/(sp_wid)
!            damp = 1.0/TNUDG * max(0.0,dummy1)
! 
!            else ; damp=0.0
!            endif


!      else
!     if (G%geoLonT(i,j) >= 0.0 .AND. G%geoLonT(i,j) <= sp_wid) then
!            if (0<=ytil .AND. ytil<=land .AND. 0<=xtil .AND. xtil<=sp_wid) then 
!
!            ! 1 / day
!            dummy1=(xtil)/(sp_wid)
!            damp = 1.0/TNUDG * max(0.0,dummy1)
!
!            elseif (0<=ytil .AND. ytil<=land .AND. (shw-sp_wid)<=xtil .AND. xtil<=shw) then
!            ! 1 / day
!            dummy1=(xtil-(shw-sp_wid))/(sp_wid)
!            damp = 1.0/TNUDG * max(0.0,dummy1)
!           
!            else ; damp=0.0
!            endif
!
            if (0<=ytil .AND. ytil<=sp_wid) then

            ! 1 / day
            dummy1=(ytil)/(sp_wid)
            damp = 1.0/TNUDG * max(0.0,dummy1)

            endif

            if (0<=xtil .AND. xtil<=sp_wid) then 
!
!            ! 1 / day
            dummy1=(xtil)/(sp_wid)
            damp = 1.0/TNUDG * max(0.0,dummy1)

            elseif ((shw-sp_wid)<=xtil .AND. xtil<=shw) then
            ! 1 / day
            dummy1=(xtil-(shw-sp_wid))/(sp_wid)
            damp = 1.0/TNUDG * max(0.0,dummy1)
           
            else ; damp=0.0
            endif 


  ! convert to 1 / seconds
      if (G%bathyT(i,j) > min_depth) then
          Idamp(i,j) = damp/86400.0
      else ; Idamp(i,j) = 0.0 ; endif

!      endif
   enddo ; enddo

  ! Compute min/max density using T_SUR/S_SUR and T_BOT/S_BOT
  call calculate_density(t_sur,s_sur,0.0,rho_sur,tv%eqn_of_state)
  !write (*,*)'Surface density in sponge:', rho_sur
  call calculate_density(t_bot,s_bot,0.0,rho_bot,tv%eqn_of_state)
  !write (*,*)'Bottom density in sponge:', rho_bot
  rho_range = rho_bot - rho_sur
  !write (*,*)'Density range in sponge:', rho_range

  if (use_ALE) then

    select case ( coordinateMode(verticalCoordinate) )

!     case ( REGRIDDING_RHO )

     case ( REGRIDDING_ZSTAR, REGRIDDING_SIGMA_SHELF_ZSTAR )   ! Initial thicknesses for z coordinates
!      do j=js,je ; do i=is,ie
!        eta1D(nz+1) = -1.0*G%bathyT(i,j)
!        do k=nz,1,-1
!          eta1D(k) =  -G%max_depth * real(k-1) / real(nz)
!          if (eta1D(k) < (eta1D(k+1) + min_thickness)) then
!            eta1D(k) = eta1D(k+1) + min_thickness
!            h(i,j,k) = min_thickness
!          else
!            h(i,j,k) = eta1D(k) - eta1D(k+1)
!          endif
!         enddo
!      enddo ; enddo

!    S_range = S_range / G%max_depth ! Convert S_range into dS/dz
!    T_range = T_range / G%max_depth ! Convert T_range into dT/dz
    do j=js,je ; do i=is,ie
       z0=0.0
       do k=1,nz
           z = z0 + 0.5*h(i,j,k) * GV%H_to_m
           if (z<zc-buff) then
               T(i,j,k) = t1
               S(i,j,k) = s1
           elseif (z>=zc-buff .and. z<=zc+buff) then
               T(i,j,k) = (t2-t1)/(2*buff)*(z-(zc-buff))+t1
               S(i,j,k) = (s2-s1)/(2*buff)*(z-(zc-buff))+s1
           elseif (z>zc+buff) then
               T(i,j,k) = t2
               S(i,j,k) = s2
           endif
           z0 = z0 +  h(i,j,k) * GV%H_to_m
        enddo
        z = z0 + 0.5*h(i,j,k) * GV%H_to_m
    enddo; enddo

!      case ( REGRIDDING_SIGMA )             ! Initial thicknesses for sigma coordinates

      case default
         call MOM_error(FATAL,"Fjord_initialize_sponges: "// &
         "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

    end select
    !  This call sets up the damping rates and interface heights.
    !  This sets the inverse damping timescale fields in the sponges.
    call initialize_ALE_sponge(Idamp, G, PF, ACSp, h, nz)

!    S_range = S_range / G%max_depth ! Convert S_range into dS/dz
!    T_range = T_range / G%max_depth ! Convert T_range into dT/dz
!    do j=js,je ; do i=is,ie
!       z0=0.0
!       do k=1,nz
!           z = z0 + 0.5*h(i,j,k) * GV%H_to_m
!           if (z<zc-buff) then
!               T(i,j,k) = t1
!               S(i,j,k) = s1
!           elseif (z>=zc-buff .and. z<=zc+buff) then
!               T(i,j,k) = (t2-t1)/(2*buff)*(z-(zc-buff))+t1
!               S(i,j,k) = (s2-s1)/(2*buff)*(z-(zc-buff))+s1
!           elseif (z>zc+buff) then
!               T(i,j,k) = t2
!               S(i,j,k) = s2
!           endif
!           z0 = z0 +  h(i,j,k) * GV%H_to_m
!        enddo
!        z = z0 + 0.5*h(i,j,k) * GV%H_to_m

!        do k=1,nz
!           z = z0 + 0.5*h(i,j,k) * GV%H_to_m
!           if (z<zc-buff) T(i,j,k) = t1
!           if (z>=zc-buff .and. z<=zc+buff) T(i,j,k) = (t2-t1)/(2*buff)*(z-(zc-buff))+t1
!           if (z>zc+buff) T(i,j,k) = t2
!           S(i,j,k) = fs_sur - exp((fd - z)/fs_L)/fs_sc + fs_m*z
!           z0 = z0 +  h(i,j,k) * GV%H_to_m
!        enddo
!   enddo ; enddo
    ! for debugging
    !i=G%iec; j=G%jec
    !do k = 1,nz
    !  call calculate_density(T(i,j,k),S(i,j,k),0.0,rho_tmp,tv%eqn_of_state)
    !  write(*,*) 'Sponge - k,h,T,S,rho,Rlay',k,h(i,j,k),T(i,j,k),S(i,j,k),rho_tmp,GV%Rlay(k)
    !enddo

    !   Now register all of the fields which are damped in the sponge.   !
    ! By default, momentum is advected vertically within the sponge, but !
    ! momentum is typically not damped within the sponge.                !

    !  The remaining calls to set_up_sponge_field can be in any order. !
    if ( associated(tv%T) ) then
      call set_up_ALE_sponge_field(T,G,tv%T,ACSp)
    endif
    if ( associated(tv%S) ) then
      call set_up_ALE_sponge_field(S,G,tv%S,ACSp)
    endif

  else ! layer mode
     do j=js,je ; do i=is,ie
       z0=0.0
       do k=1,nz
           z = z0 + 0.5*h(i,j,k) * GV%H_to_m
           if (z<zc-buff) then
               T(i,j,k) = t1
               S(i,j,k) = s1
           elseif (z>=zc-buff .and. z<=zc+buff) then
               T(i,j,k) = (t2-t1)/(2*buff)*(z-(zc-buff))+t1
               S(i,j,k) = (s2-s1)/(2*buff)*(z-(zc-buff))+s1
           elseif (z>zc+buff) then
               T(i,j,k) = t2
               S(i,j,k) = s2
           endif
           z0 = z0 +  h(i,j,k) * GV%H_to_m
        enddo
        z = z0 + 0.5*h(i,j,k) * GV%H_to_m
    enddo; enddo

    
!       ! 1) Read eta, salt and temp from IC file
!       call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
!       inputdir = slasher(inputdir)
       ! GM: get two different files, one with temp and one with salt values
       ! this is work around to avoid having wrong values near the surface
       ! because of the FIT_SALINITY option. To get salt values right in the
       ! sponge, FIT_SALINITY=False. The oposite is true for temp. One can
       ! combined the *correct* temp and salt values in one file instead.
!       call get_param(PF, mdl, "Fjord_SPONGE_FILE", state_file, &
!                 "The name of the file with temps., salts. and interfaces to \n"// &
!                 " damp toward.", fail_if_missing=.true.)
!       call get_param(PF, mdl, "SPONGE_PTEMP_VAR", temp_var, &
!                 "The name of the potential temperature variable in \n"//&
!                 "SPONGE_STATE_FILE.", default="Temp")
!       call get_param(PF, mdl, "SPONGE_SALT_VAR", salt_var, &
!                 "The name of the salinity variable in \n"//&
!                 "SPONGE_STATE_FILE.", default="Salt")
!       call get_param(PF, mdl, "SPONGE_ETA_VAR", eta_var, &
!                 "The name of the interface height variable in \n"//&
!                 "SPONGE_STATE_FILE.", default="eta")
!
!       !read temp and eta
!       filename = trim(inputdir)//trim(state_file)
!       if (.not.file_exists(filename, G%Domain)) &
!          call MOM_error(FATAL, " Fjord_initialize_sponges: Unable to open "//trim(filename))
!       call MOM_read_data(filename, eta_var, eta(:,:,:), G%Domain)
!       call MOM_read_data(filename, temp_var, T(:,:,:), G%Domain)
!       call MOM_read_data(filename, salt_var, S(:,:,:), G%Domain)
!
!       ! for debugging
!       !i=G%iec; j=G%jec
!       !do k = 1,nz
!       !    call calculate_density(T(i,j,k),S(i,j,k),0.0,rho_tmp,tv%eqn_of_state)
!       !    write(*,*) 'Sponge - k,eta,T,S,rho,Rlay',k,eta(i,j,k),T(i,j,k),&
       !                S(i,j,k),rho_tmp,GV%Rlay(k)
       !enddo

       ! Set the inverse damping rates so that the model will know where to
       ! apply the sponges, along with the interface heights.
       call initialize_sponge(Idamp, eta, G, PF, CSp)
!       ! Apply sponge in tracer fields
       call set_up_sponge_field(T, tv%T, G, nz, CSp)
       call set_up_sponge_field(S, tv%S, G, nz, CSp)

  endif

end subroutine Fjord_initialize_sponges

!> This subroutine sets the properties of flow at open boundary conditions.
!! This is copied and modified from the DOME inflow example. 
subroutine Fjord_set_OBC_data(OBC, tv, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv  !< A structure containing pointers to any
                              !! available thermodynamic fields, including potential
                              !! temperature and salinity or mixed layer density. Absent
                              !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                              !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

! Local variables
  ! The following variables are used to set the target temperature and salinity.
  real :: T0(SZK_(G)), S0(SZK_(G))
  real :: pres(SZK_(G))      ! An array of the reference pressure in Pa.
  real :: drho_dT(SZK_(G))   ! Derivative of density with temperature in kg m-3 K-1.                              !
  real :: drho_dS(SZK_(G))   ! Derivative of density with salinity in kg m-3 PSU-1.                             !
  real :: rho_guess(SZK_(G)) ! Potential density at T0 & S0 in kg m-3.
  ! The following variables are used to set up the transport
  real :: tr_0, y1, y2, tr_k, rst, rsb, rc, v_k, lon_im1
  real :: D_edge            ! The thickness in m of the dense fluid at the
                            ! inner edge of the inflow.
  real :: g_prime_tot       ! The reduced gravity across all layers, m s-2.
  real :: Def_Rad           ! The deformation radius, based on fluid of
                            ! thickness D_edge, in the same units as lat.
  real :: Ri_trans          ! The shear Richardson number in the transition
                            ! region of the specified shear profile.
  character(len=40)  :: mdl = "Fjord_set_OBC_data" ! This subroutine's name.
  character(len=32)  :: name
  integer :: i, j, k, itt, is, ie, js, je, isd, ied, jsd, jed, m, nz, NTR
  integer :: IsdB, IedB, JsdB, JedB
  type(OBC_segment_type), pointer :: segment
  type(tracer_type), pointer      :: tr_ptr

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! The following variables should be transformed into runtime parameters.
! D_edge = 10.0  ! The thickness of dense fluid in the inflow.
! Ri_trans = 1.0/3.0 ! The shear Richardson number in the transition region
                     ! region of the specified shear profile.

  if (.not.associated(OBC)) return

! g_prime_tot = (GV%g_Earth/GV%Rho0)*2.0
! Def_Rad = sqrt(D_edge*g_prime_tot) / (1.0e-4*1000.0)
! tr_0 = (-D_edge*sqrt(D_edge*g_prime_tot)*0.5e3*Def_Rad) * GV%m_to_H

  if (OBC%number_of_segments /= 1) then
    print *, 'Error in Fjord OBC segment setup'
    return   !!! Need a better error message here
  endif
  segment => OBC%segment(1)
  if (.not. segment%on_pe) return

  NTR = tr_Reg%NTR !Number of tracers?
  allocate(segment%field(NTR)) 

! Here add select case for river OBC or just shelf forcing OBC

!  do k=1
   k=1
   ! New way
    isd = segment%HI%isd ; ied = segment%HI%ied
    JsdB = segment%HI%JsdB ; JedB = segment%HI%JedB
    do J=JsdB,JedB ; do i=isd,ied
      segment%normal_trans(i,J,k) = 0.01;
      segment%normal_vel(i,J,k) = 0.01;
    enddo ; enddo
!  enddo

  !   The inflow values of temperature and salinity also need to be set here if
  ! these variables are used.  The following code is just a naive example.
  if (associated(tv%S)) then
    ! In this example, all S inflows have values of 35 psu.
    name = 'salt'
    call tracer_name_lookup(tr_Reg, tr_ptr, name)
    call register_segment_tracer(tr_ptr, param_file, GV, segment, OBC_scalar=0.0)
  endif
  if (associated(tv%T)) then
    name = 'temp'
    call tracer_name_lookup(tr_Reg, tr_ptr, name)
    call register_segment_tracer(tr_ptr, param_file, GV, segment, OBC_scalar=0.0)
  endif

  ! Dye tracers - fight with T,S???
  ! First dye - only one with OBC values
  ! This field(1) requires tr_D1 to be the first tracer.
!  allocate(segment%field(1)%buffer_src(segment%HI%isd:segment%HI%ied,segment%HI%JsdB:segment%HI%JedB,nz))
!  do k=1,nz ; do j=segment%HI%jsd,segment%HI%jed ; do i=segment%HI%isd,segment%HI%ied
!    if (k < nz/2) then ; segment%field(1)%buffer_src(i,j,k) = 0.0
!    else ; segment%field(1)%buffer_src(i,j,k) = 1.0 ; endif
!  enddo ; enddo ; enddo
!  name = 'tr_D1'
!  call tracer_name_lookup(tr_Reg, tr_ptr, name)
!  call register_segment_tracer(tr_ptr, param_file, GV, &
!                               OBC%segment(1), OBC_array=.true.)

  ! All tracers but the first have 0 concentration in their inflows. As this
  ! is the default value, the following calls are unnecessary.
!  do m=2,NTR
!    if (m < 10) then ; write(name,'("tr_D",I1.1)') m
!    else ; write(name,'("tr_D",I2.2)') m ; endif
!    call tracer_name_lookup(tr_Reg, tr_ptr, name)
!    call register_segment_tracer(tr_ptr, param_file, GV, &
!                                 OBC%segment(1), OBC_scalar=0.0)
!  enddo

end subroutine Fjord_set_OBC_data

!> \namespace isomip_initialization
!!
!!  The module configures the ISOMIP test case.
end module Fjord_initialization

