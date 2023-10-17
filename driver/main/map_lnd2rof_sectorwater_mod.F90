module map_lnd2rof_sectorwater_mod

    !---------------------------------------------------------------------
    !
    ! Purpose:
    !
    ! This module contains routines for mapping the sectoral water fields from the LND grid onto
    ! the ROF grid.
    !
    ! These routines could go in prep_rof_mod, but are separated into their own module for
    ! the sake of (1) testability: this module has fewer dependencies than prep_rof_mod;
    ! and (2) symmetry with the lnd2glc and glc2lnd custom mapping routines, which also
    ! have their own modules.
  
  #include "shr_assert.h"
    use shr_kind_mod, only : r8 => shr_kind_r8
    use mct_mod
    use seq_map_type_mod, only : seq_map
    use seq_map_mod, only : seq_map_map
  
    implicit none
    private
  
    ! ------------------------------------------------------------------------
    ! Public interfaces
    ! ------------------------------------------------------------------------
  
    public :: map_lnd2rof_sectorwater  ! map sector water fluxes from lnd -> rof grid
  
    ! ------------------------------------------------------------------------
    ! Private interfaces
    ! ------------------------------------------------------------------------
  
    private :: map_rof2lnd_volr  ! map volr from rof -> lnd grid
  
    character(len=*), parameter, private :: sourcefile = &
         __FILE__
  
  contains
  
    subroutine map_lnd2rof_sectorwater(l2r_l, r2x_r, dom_withd_flux_field, dom_rf_flux_field, &
         liv_withd_flux_field, liv_rf_flux_field, elec_withd_flux_field, elec_rf_flux_field,  &
         mfc_withd_flux_field, mfc_rf_flux_field, min_withd_flux_field, min_rf_flux_field,    &
         avwts_s, avwtsfld_s, mapper_Fl2r, mapper_Fr2l, l2r_r)
      !---------------------------------------------------------------
      ! Description
      ! Do custom mapping for the sectoral water fluxex, from land -> rof.
      !
      ! The basic idea is that we want to pull/add sectoral water fluxes out of/in ROF cells proportionally to
      ! the river volume (volr) in each cell. This is important in cases where the various
      ! ROF cells overlapping a CLM cell have very different volr: If we didn't do this
      ! volr-normalized remapping, we'd try to extract the same amount of water from each
      ! of the ROF cells, which would be more likely to have withdrawals exceeding
      ! available volr.
      !
      ! (Both RTM and MOSART have code to handle excess withdrawals, by pulling the excess
      ! directly out of the ocean, but we'd like to avoid resorting to this as much as
      ! possible.)
      !
      ! This mapping works by:
      !
      ! (1) Normalizing the land's sector water flux by volr
      !
      ! (2) Mapping this volr-normalized flux to the rof grid
      !
      ! (3) Converting the mapped, volr-normalized flux back to a normal
      !     (non-volr-normalized) flux on the rof grid.
      !
      ! This assumes that the following fields are contained in the attribute vector
      ! arguments:
      !
      ! - l2r_l: field given by sectorX_withd/rf_flux_field (read)
      ! - l2r_r: field given by sectorX_withd/rf_flux_field (set)
      ! - r2x_r: 'Flrr_volrmch' (read)
      !
      ! Arguments
      type(mct_aVect)  , intent(in)    :: l2r_l            ! lnd -> rof fields on the land grid
      type(mct_aVect)  , intent(in)    :: r2x_r            ! rof -> cpl fields on the rof grid
      character(len=*) , intent(in)    :: dom_withd_flux_field ! name of domestic withdrawal field to remap
      character(len=*) , intent(in)    :: dom_rf_flux_field ! name of domestic return field to remap
      character(len=*) , intent(in)    :: liv_withd_flux_field ! name of livestock withdrawal field to remap
      character(len=*) , intent(in)    :: liv_rf_flux_field ! name of livestock return field to remap
      character(len=*) , intent(in)    :: elec_withd_flux_field ! name of thermoelectri withdrawal field to remap
      character(len=*) , intent(in)    :: elec_rf_flux_field ! name of thermoelectri return field to remap
      character(len=*) , intent(in)    :: mfc_withd_flux_field ! name of manufacturing withdrawal field to remap
      character(len=*) , intent(in)    :: mfc_rf_flux_field ! name of manufacturing return field to remap
      character(len=*) , intent(in)    :: min_withd_flux_field ! name of mining withdrawal field to remap
      character(len=*) , intent(in)    :: min_rf_flux_field ! name of mining return field to remap
      
  
      type(mct_aVect)  , intent(in)    :: avwts_s          ! attr vect for source weighting
      character(len=*) , intent(in)    :: avwtsfld_s       ! field in avwts_s to use
      type(seq_map)    , intent(inout) :: mapper_Fl2r      ! flux mapper for mapping lnd -> rof
      type(seq_map)    , intent(inout) :: mapper_Fr2l      ! flux mapper for mapping rof -> lnd
      type(mct_aVect)  , intent(inout) :: l2r_r            ! lnd -> rof fields on the rof grid
      !
      ! Local variables
      integer :: r, l
      integer :: lsize_l  ! number of land points
      integer :: lsize_r  ! number of rof points
      type(mct_avect) :: dom_withd_l_av  ! temporary attribute vector holding domestic withdrawal fluxes on the land grid
      type(mct_avect) :: dom_withd_r_av  ! temporary attribute vector holding domestic withdrawal fluxes on the rof grid
      type(mct_avect) :: dom_rf_l_av  ! temporary attribute vector holding domestic return flow fluxes on the land grid
      type(mct_avect) :: dom_rf_r_av  ! temporary attribute vector holding domestic return flow fluxes on the rof grid
      
      type(mct_avect) :: liv_withd_l_av  ! temporary attribute vector holding livestock withdrawal fluxes on the land grid
      type(mct_avect) :: liv_withd_r_av  ! temporary attribute vector holding livestock withdrawal fluxes on the rof grid
      type(mct_avect) :: liv_rf_l_av  ! temporary attribute vector holding livestock return flow fluxes on the land grid
      type(mct_avect) :: liv_rf_r_av  ! temporary attribute vector holding livestock return flow fluxes on the rof grid
  
      type(mct_avect) :: elec_withd_l_av  ! temporary attribute vector holding thermoelectric withdrawal fluxes on the land grid
      type(mct_avect) :: elec_withd_r_av  ! temporary attribute vector holding thermoelectric withdrawal fluxes on the rof grid
      type(mct_avect) :: elec_rf_l_av  ! temporary attribute vector holding thermoelectric return flow fluxes on the land grid
      type(mct_avect) :: elec_rf_r_av  ! temporary attribute vector holding thermoelectric return flow fluxes on the rof grid
  
      type(mct_avect) :: mfc_withd_l_av  ! temporary attribute vector holding manufacturing withdrawal fluxes on the land grid
      type(mct_avect) :: mfc_withd_r_av  ! temporary attribute vector holding manufacturing withdrawal fluxes on the rof grid
      type(mct_avect) :: mfc_rf_l_av  ! temporary attribute vector holding manufacturing return flow fluxes on the land grid
      type(mct_avect) :: mfc_rf_r_av  ! temporary attribute vector holding manufacturing return flow fluxes on the rof grid
  
      type(mct_avect) :: min_withd_l_av  ! temporary attribute vector holding mining withdrawal fluxes on the land grid
      type(mct_avect) :: min_withd_r_av  ! temporary attribute vector holding mining withdrawal fluxes on the rof grid
      type(mct_avect) :: min_rf_l_av  ! temporary attribute vector holding mining return flow fluxes on the land grid
      type(mct_avect) :: min_rf_r_av  ! temporary attribute vector holding mining return flow fluxes on the rof grid
  
      ! The following need to be pointers to satisfy the MCT interface:
      real(r8), pointer :: volr_r(:)             ! river volume on the rof grid
      real(r8), pointer :: volr_l(:)             ! river volume on the land grid
  
      real(r8), pointer :: dom_withd_flux_l(:)       ! domestic withdrawal flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: dom_withd_flux_r(:)       ! domestic withdrawal flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: dom_withd_normalized_l(:) ! domestic withdrawal normalized by volr, land grid
      real(r8), pointer :: dom_withd_normalized_r(:) ! domestic withdrawal normalized by volr, rof grid
      real(r8), pointer :: dom_withd_volr0_l(:)      ! domestic withdrawal where volr <= 0, land grid
      real(r8), pointer :: dom_withd_volr0_r(:)      ! domestic withdrawal where volr <= 0, rof grid
  
      real(r8), pointer :: dom_rf_flux_l(:)       ! domestic return flow flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: dom_rf_flux_r(:)       ! domestic return flow flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: dom_rf_normalized_l(:) ! domestic return flow normalized by volr, land grid
      real(r8), pointer :: dom_rf_normalized_r(:) ! domestic return flow normalized by volr, rof grid
      real(r8), pointer :: dom_rf_volr0_l(:)      ! domestic return flow where volr <= 0, land grid
      real(r8), pointer :: dom_rf_volr0_r(:)      ! domestic return flow where volr <= 0, rof grid
  
      real(r8), pointer :: liv_withd_flux_l(:)       ! livestock withdrawal flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: liv_withd_flux_r(:)       ! livestock withdrawal flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: liv_withd_normalized_l(:) ! livestock withdrawal normalized by volr, land grid
      real(r8), pointer :: liv_withd_normalized_r(:) ! livestock withdrawal normalized by volr, rof grid
      real(r8), pointer :: liv_withd_volr0_l(:)      ! livestock withdrawal where volr <= 0, land grid
      real(r8), pointer :: liv_withd_volr0_r(:)      ! livestock withdrawal where volr <= 0, rof grid
  
      real(r8), pointer :: liv_rf_flux_l(:)       ! livestock return flow flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: liv_rf_flux_r(:)       ! livestock return flow flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: liv_rf_normalized_l(:) ! livestock return flow normalized by volr, land grid
      real(r8), pointer :: liv_rf_normalized_r(:) ! livestock return flow normalized by volr, rof grid
      real(r8), pointer :: liv_rf_volr0_l(:)      ! livestock return flow where volr <= 0, land grid
      real(r8), pointer :: liv_rf_volr0_r(:)      ! livestock return flow where volr <= 0, rof grid
      
  
      real(r8), pointer :: elec_withd_flux_l(:)       ! thermoelectric withdrawal flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: elec_withd_flux_r(:)       ! thermoelectric withdrawal flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: elec_withd_normalized_l(:) ! thermoelectric withdrawal normalized by volr, land grid
      real(r8), pointer :: elec_withd_normalized_r(:) ! thermoelectric withdrawal normalized by volr, rof grid
      real(r8), pointer :: elec_withd_volr0_l(:)      ! thermoelectric withdrawal where volr <= 0, land grid
      real(r8), pointer :: elec_withd_volr0_r(:)      ! thermoelectric withdrawal where volr <= 0, rof grid
  
      real(r8), pointer :: elec_rf_flux_l(:)       ! thermoelectric return flow flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: elec_rf_flux_r(:)       ! thermoelectric return flow flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: elec_rf_normalized_l(:) ! thermoelectric return flow normalized by volr, land grid
      real(r8), pointer :: elec_rf_normalized_r(:) ! thermoelectric return flow normalized by volr, rof grid
      real(r8), pointer :: elec_rf_volr0_l(:)      ! thermoelectric return flow where volr <= 0, land grid
      real(r8), pointer :: elec_rf_volr0_r(:)      ! thermoelectric return flow where volr <= 0, rof grid
      
  
      real(r8), pointer :: mfc_withd_flux_l(:)       ! manufacturing withdrawal flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: mfc_withd_flux_r(:)       ! manufacturing withdrawal flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: mfc_withd_normalized_l(:) ! manufacturing withdrawal normalized by volr, land grid
      real(r8), pointer :: mfc_withd_normalized_r(:) ! manufacturing withdrawal normalized by volr, rof grid
      real(r8), pointer :: mfc_withd_volr0_l(:)      ! manufacturing withdrawal where volr <= 0, land grid
      real(r8), pointer :: mfc_withd_volr0_r(:)      ! manufacturing withdrawal where volr <= 0, rof grid
  
      real(r8), pointer :: mfc_rf_flux_l(:)       ! manufacturing return flow flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: mfc_rf_flux_r(:)       ! manufacturing return flow flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: mfc_rf_normalized_l(:) ! manufacturing return flow normalized by volr, land grid
      real(r8), pointer :: mfc_rf_normalized_r(:) ! manufacturing return flow normalized by volr, rof grid
      real(r8), pointer :: mfc_rf_volr0_l(:)      ! manufacturing return flow where volr <= 0, land grid
      real(r8), pointer :: mfc_rf_volr0_r(:)      ! manufacturing return flow where volr <= 0, rof grid
      
  
      real(r8), pointer :: min_withd_flux_l(:)       ! mining withdrawal flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: min_withd_flux_r(:)       ! mining withdrawal flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: min_withd_normalized_l(:) ! mining withdrawal normalized by volr, land grid
      real(r8), pointer :: min_withd_normalized_r(:) ! mining withdrawal normalized by volr, rof grid
      real(r8), pointer :: min_withd_volr0_l(:)      ! mining withdrawal where volr <= 0, land grid
      real(r8), pointer :: min_withd_volr0_r(:)      ! mining withdrawal where volr <= 0, rof grid
  
      real(r8), pointer :: min_rf_flux_l(:)       ! mining return flow flux on the land grid [kg m-2 s-1]
      real(r8), pointer :: min_rf_flux_r(:)       ! mining return flow flux on the rof grid [kg m-2 s-1]
      real(r8), pointer :: min_rf_normalized_l(:) ! mining return flow normalized by volr, land grid
      real(r8), pointer :: min_rf_normalized_r(:) ! mining return flow normalized by volr, rof grid
      real(r8), pointer :: min_rf_volr0_l(:)      ! mining return flow where volr <= 0, land grid
      real(r8), pointer :: min_rf_volr0_r(:)      ! mining return flow where volr <= 0, rof grid
      
      character(len=*), parameter :: volr_field             = 'Flrr_volrmch'
  
      character(len=*), parameter :: dom_withd_normalized_field = 'Flrl_dom_withd_normalized'
      character(len=*), parameter :: dom_withd_volr0_field      = 'Flrl_dom_withd_volr0'
      character(len=*), parameter :: fields_to_remap_dom_withd = dom_withd_normalized_field // ':' // dom_withd_volr0_field
  
      character(len=*), parameter :: dom_rf_normalized_field = 'Flrl_dom_rf_normalized'
      character(len=*), parameter :: dom_rf_volr0_field      = 'Flrl_dom_rf_volr0'
      character(len=*), parameter :: fields_to_remap_dom_rf = dom_rf_normalized_field // ':' // dom_rf_volr0_field
  
      character(len=*), parameter :: liv_withd_normalized_field = 'Flrl_liv_withd_normalized'
      character(len=*), parameter :: liv_withd_volr0_field      = 'Flrl_liv_withd_volr0'
      character(len=*), parameter :: fields_to_remap_liv_withd = liv_withd_normalized_field // ':' // liv_withd_volr0_field
  
      character(len=*), parameter :: liv_rf_normalized_field = 'Flrl_liv_rf_normalized'
      character(len=*), parameter :: liv_rf_volr0_field      = 'Flrl_liv_rf_volr0'
      character(len=*), parameter :: fields_to_remap_lif_rf = liv_rf_normalized_field // ':' // liv_rf_volr0_field
  
      character(len=*), parameter :: elec_withd_normalized_field = 'Flrl_elec_withd_normalized'
      character(len=*), parameter :: elec_withd_volr0_field      = 'Flrl_elec_withd_volr0'
      character(len=*), parameter :: fields_to_remap_elec_withd = elec_withd_normalized_field // ':' // elec_withd_volr0_field
  
      character(len=*), parameter :: elec_rf_normalized_field = 'Flrl_elec_rf_normalized'
      character(len=*), parameter :: elec_rf_volr0_field      = 'Flrl_elec_rf_volr0'
      character(len=*), parameter :: fields_to_remap_elec_rf = elec_rf_normalized_field // ':' // elec_rf_volr0_field
  
      character(len=*), parameter :: mfc_withd_normalized_field = 'Flrl_mfc_withd_normalized'
      character(len=*), parameter :: mfc_withd_volr0_field      = 'Flrl_mfc_withd_volr0'
      character(len=*), parameter :: fields_to_remap_mfc_withd = mfc_withd_normalized_field // ':' // mfc_withd_volr0_field
  
      character(len=*), parameter :: mfc_rf_normalized_field = 'Flrl_mfc_rf_normalized'
      character(len=*), parameter :: mfc_rf_volr0_field      = 'Flrl_mfc_rf_volr0'
      character(len=*), parameter :: fields_to_remap_mfc_rf = mfc_rf_normalized_field // ':' // mfc_rf_volr0_field
  
      character(len=*), parameter :: min_withd_normalized_field = 'Flrl_min_withd_normalized'
      character(len=*), parameter :: min_withd_volr0_field      = 'Flrl_min_withd_volr0'
      character(len=*), parameter :: fields_to_remap_min_withd = min_withd_normalized_field // ':' // min_withd_volr0_field
  
      character(len=*), parameter :: min_rf_normalized_field = 'Flrl_min_rf_normalized'
      character(len=*), parameter :: min_rf_volr0_field      = 'Flrl_min_rf_volr0'
      character(len=*), parameter :: fields_to_remap_min_rf = min_rf_normalized_field // ':' // min_rf_volr0_field
  
      
      !---------------------------------------------------------------
  
      ! ------------------------------------------------------------------------
      ! Determine attribute vector sizes
      ! ------------------------------------------------------------------------
  
      lsize_l = mct_aVect_lsize(l2r_l)
      lsize_r = mct_aVect_lsize(l2r_r)
  
      ! ------------------------------------------------------------------------
      ! Extract the necessary fields from attribute vectors
      ! ------------------------------------------------------------------------
  
      allocate(dom_withd_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, dom_withd_flux_field, dom_withd_flux_l)
      allocate(dom_rf_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, dom_rf_flux_field, dom_rf_flux_l)
  
      allocate(liv_withd_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, liv_withd_flux_field, liv_withd_flux_l)
      allocate(liv_rf_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, liv_rf_flux_field, liv_rf_flux_l)
  
      allocate(elec_withd_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, elec_withd_flux_field, elec_withd_flux_l)
      allocate(elec_rf_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, elec_rf_flux_field, elec_rf_flux_l)
  
      allocate(mfc_withd_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, mfc_withd_flux_field, mfc_withd_flux_l)
      allocate(mfc_rf_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, mfc_rf_flux_field, mfc_rf_flux_l)
  
      allocate(min_withd_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, min_withd_flux_field, min_withd_flux_l)
      allocate(min_rf_flux_l(lsize_l))
      call mct_aVect_exportRattr(l2r_l, min_rf_flux_field, min_rf_flux_l)
  
      allocate(volr_r(lsize_r))
      call mct_aVect_exportRattr(r2x_r, volr_field, volr_r)
  
      ! ------------------------------------------------------------------------
      ! Adjust volr_r, and map it to the land grid
      ! ------------------------------------------------------------------------
  
      ! Treat any rof point with volr < 0 as if it had volr = 0. Negative volr values can
      ! arise in RTM. This fix is needed to avoid mapping negative sector water fluxes to those
      ! cells: while conservative, this would be unphysical (it would mean that sectoral water fluxes
      ! actually adds water to those cells).
      do r = 1, lsize_r
         if (volr_r(r) < 0._r8) then
            volr_r(r) = 0._r8
         end if
      end do
  
      allocate(volr_l(lsize_l))
      call map_rof2lnd_volr(volr_r, mapper_Fr2l, volr_l)
  
      ! ------------------------------------------------------------------------
      ! Determine sector water fluxes normalized by volr
      !
      ! In order to avoid possible divide by 0, as well as to handle non-sensical negative
      ! volr on the land grid, we divide the land's sectoral water flux into two separate flux
      ! components: a component where we have positive volr on the land grid (put in
      ! sectorX_withd/rf_normalized_l, which is mapped using volr-normalization) and a component where
      ! we have zero or negative volr on the land grid (put in sectorX_wihtd/rf_volr0_l, which is
      ! mapped as a standard flux). We then remap both of these components to the rof grid,
      ! and then finally add the two components to determine the total sector water flux on
      ! the rof grid.
      ! ------------------------------------------------------------------------
  
      allocate(dom_withd_normalized_l(lsize_l))
      allocate(dom_withd_volr0_l(lsize_l))
      allocate(dom_rf_normalized_l(lsize_l))
      allocate(dom_rf_volr0_l(lsize_l))
  
      allocate(liv_withd_normalized_l(lsize_l))
      allocate(liv_withd_volr0_l(lsize_l))
      allocate(liv_rf_normalized_l(lsize_l))
      allocate(liv_rf_volr0_l(lsize_l))
  
      allocate(elec_withd_normalized_l(lsize_l))
      allocate(elec_withd_volr0_l(lsize_l))
      allocate(elec_rf_normalized_l(lsize_l))
      allocate(elec_rf_volr0_l(lsize_l))
  
      allocate(mfc_withd_normalized_l(lsize_l))
      allocate(mfc_withd_volr0_l(lsize_l))
      allocate(mfc_rf_normalized_l(lsize_l))
      allocate(mfc_rf_volr0_l(lsize_l))
  
      allocate(min_withd_normalized_l(lsize_l))
      allocate(min_withd_volr0_l(lsize_l))
      allocate(min_rf_normalized_l(lsize_l))
      allocate(min_rf_volr0_l(lsize_l))
  
      do l = 1, lsize_l
         if (volr_l(l) > 0._r8) then
            dom_withd_normalized_l(l) = dom_withd_flux_l(l) / volr_l(l)
            dom_withd_volr0_l(l)      = 0._r8
            liv_withd_normalized_l(l) = liv_withd_flux_l(l) / volr_l(l)
            liv_withd_volr0_l(l)      = 0._r8
            elec_withd_normalized_l(l) = elec_withd_flux_l(l) / volr_l(l)
            elec_withd_volr0_l(l)      = 0._r8
            mfc_withd_normalized_l(l) = mfc_withd_flux_l(l) / volr_l(l)
            mfc_withd_volr0_l(l)      = 0._r8
            min_withd_normalized_l(l) = min_withd_flux_l(l) / volr_l(l)
            min_withd_volr0_l(l)      = 0._r8
  
            dom_rf_normalized_l(l) = dom_rf_flux_l(l) / volr_l(l)
            dom_rf_volr0_l(l)      = 0._r8
            liv_rf_normalized_l(l) = liv_rf_flux_l(l) / volr_l(l)
            liv_rf_volr0_l(l)      = 0._r8
            elec_rf_normalized_l(l) = elec_rf_flux_l(l) / volr_l(l)
            elec_rf_volr0_l(l)      = 0._r8
            mfc_rf_normalized_l(l) = mfc_rf_flux_l(l) / volr_l(l)
            mfc_rf_volr0_l(l)      = 0._r8
            min_rf_normalized_l(l) = min_rf_flux_l(l) / volr_l(l)
            min_rf_volr0_l(l)      = 0._r8
         else
            dom_withd_normalized_l(l) = 0._r8
            dom_withd_volr0_l(l)      = dom_withd_flux_l(l)
            liv_withd_normalized_l(l) = 0._r8
            liv_withd_volr0_l(l)      = liv_withd_flux_l(l)
            elec_withd_normalized_l(l) = 0._r8
            elec_withd_volr0_l(l)      = elec_withd_flux_l(l)
            mfc_withd_normalized_l(l) = 0._r8
            mfc_withd_volr0_l(l)      = mfc_withd_flux_l(l)
            min_withd_normalized_l(l) = 0._r8
            min_withd_volr0_l(l)      = min_withd_flux_l(l)
        
            dom_rf_normalized_l(l) = 0._r8
            dom_rf_volr0_l(l)      = dom_rf_flux_l(l)
            liv_rf_normalized_l(l) = 0._r8
            liv_rf_volr0_l(l)      = liv_rf_flux_l(l)
            elec_rf_normalized_l(l) = 0._r8
            elec_rf_volr0_l(l)      = elec_rf_flux_l(l)
            mfc_rf_normalized_l(l) = 0._r8
            mfc_rf_volr0_l(l)      = mfc_rf_flux_l(l)
            min_rf_normalized_l(l) = 0._r8
            min_rf_volr0_l(l)      = min_rf_flux_l(l)
         end if
      end do
  
      ! ------------------------------------------------------------------------
      ! Map sector water fluxes
      ! ------------------------------------------------------------------------
  
      call mct_aVect_init(dom_withd_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(dom_withd_l_av, dom_withd_normalized_field, dom_withd_normalized_l)
      call mct_aVect_importRattr(dom_withd_l_av, dom_withd_volr0_field, dom_withd_volr0_l)
      call mct_aVect_init(dom_withd_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(dom_rf_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(dom_rf_l_av, dom_rf_normalized_field, dom_rf_normalized_l)
      call mct_aVect_importRattr(dom_rf_l_av, dom_rf_volr0_field, dom_rf_volr0_l)
      call mct_aVect_init(dom_rf_r_av, rList = fields_to_remap, lsize = lsize_r)
      
      call mct_aVect_init(liv_withd_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(liv_withd_l_av, liv_withd_normalized_field, liv_withd_normalized_l)
      call mct_aVect_importRattr(liv_withd_l_av, liv_withd_volr0_field, liv_withd_volr0_l)
      call mct_aVect_init(liv_withd_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(liv_rf_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(liv_rf_l_av, liv_rf_normalized_field, liv_rf_normalized_l)
      call mct_aVect_importRattr(liv_rf_l_av, liv_rf_volr0_field, liv_rf_volr0_l)
      call mct_aVect_init(liv_rf_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(elec_withd_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(elec_withd_l_av, elec_withd_normalized_field, elec_withd_normalized_l)
      call mct_aVect_importRattr(elec_withd_l_av, elec_withd_volr0_field, elec_withd_volr0_l)
      call mct_aVect_init(elec_withd_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(elec_rf_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(elec_rf_l_av, elec_rf_normalized_field, elec_rf_normalized_l)
      call mct_aVect_importRattr(elec_rf_l_av, elec_rf_volr0_field, elec_rf_volr0_l)
      call mct_aVect_init(elec_rf_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(mfc_withd_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(mfc_withd_l_av, mfc_withd_normalized_field, mfc_withd_normalized_l)
      call mct_aVect_importRattr(mfc_withd_l_av, mfc_withd_volr0_field, mfc_withd_volr0_l)
      call mct_aVect_init(mfc_withd_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(mfc_rf_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(mfc_rf_l_av, mfc_rf_normalized_field, mfc_rf_normalized_l)
      call mct_aVect_importRattr(mfc_rf_l_av, mfc_rf_volr0_field, mfc_rf_volr0_l)
      call mct_aVect_init(mfc_rf_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(min_withd_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(min_withd_l_av, min_withd_normalized_field, min_withd_normalized_l)
      call mct_aVect_importRattr(min_withd_l_av, min_withd_volr0_field, min_withd_volr0_l)
      call mct_aVect_init(min_withd_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      call mct_aVect_init(min_rf_l_av, rList = fields_to_remap, lsize = lsize_l)
      call mct_aVect_importRattr(min_rf_l_av, min_rf_normalized_field, min_rf_normalized_l)
      call mct_aVect_importRattr(min_rf_l_av, min_rf_volr0_field, min_rf_volr0_l)
      call mct_aVect_init(min_rf_r_av, rList = fields_to_remap, lsize = lsize_r)
  
      ! This mapping uses the same options (such as avwts) as is used for mapping all other
      ! fields in prep_rof_calc_l2r_rx
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = dom_withd_l_av, &
           av_d = dom_withd_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = dom_rf_l_av, &
           av_d = dom_rf_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
           
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = liv_withd_l_av, &
           av_d = liv_withd_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = liv_rf_l_av, &
           av_d = liv_rf_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
        
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = elec_withd_l_av, &
           av_d = elec_withd_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = elec_rf_l_av, &
           av_d = elec_rf_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = mfc_withd_l_av, &
           av_d = mfc_withd_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = mfc_rf_l_av, &
           av_d = mfc_rf_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = min_withd_l_av, &
           av_d = min_withd_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      call seq_map_map(mapper = mapper_Fl2r, &
           av_s = min_rf_l_av, &
           av_d = min_rf_r_av, &
           fldlist = fields_to_remap, &
           norm = .true., &
           avwts_s = avwts_s, &
           avwtsfld_s = avwtsfld_s)
  
      allocate(dom_withd_normalized_r(lsize_r))
      allocate(dom_withd_volr0_r(lsize_r))
      call mct_aVect_exportRattr(dom_withd_r_av, dom_withd_normalized_field, dom_withd_normalized_r)
      call mct_aVect_exportRattr(dom_withd_r_av, dom_withd_volr0_field, dom_withd_volr0_r)
  
      allocate(dom_rf_normalized_r(lsize_r))
      allocate(dom_rf_volr0_r(lsize_r))
      call mct_aVect_exportRattr(dom_rf_r_av, dom_rf_normalized_field, dom_rf_normalized_r)
      call mct_aVect_exportRattr(dom_rf_r_av, dom_rf_volr0_field, dom_rf_volr0_r)
  
      allocate(liv_withd_normalized_r(lsize_r))
      allocate(liv_withd_volr0_r(lsize_r))
      call mct_aVect_exportRattr(liv_withd_r_av, liv_withd_normalized_field, liv_withd_normalized_r)
      call mct_aVect_exportRattr(liv_withd_r_av, liv_withd_volr0_field, liv_withd_volr0_r)
  
      allocate(liv_rf_normalized_r(lsize_r))
      allocate(liv_rf_volr0_r(lsize_r))
      call mct_aVect_exportRattr(liv_rf_r_av, liv_rf_normalized_field, liv_rf_normalized_r)
      call mct_aVect_exportRattr(liv_rf_r_av, liv_rf_volr0_field, liv_rf_volr0_r)
  
      allocate(elec_withd_normalized_r(lsize_r))
      allocate(elec_withd_volr0_r(lsize_r))
      call mct_aVect_exportRattr(elec_withd_r_av, elec_withd_normalized_field, elec_withd_normalized_r)
      call mct_aVect_exportRattr(elec_withd_r_av, elec_withd_volr0_field, elec_withd_volr0_r)
  
      allocate(elec_rf_normalized_r(lsize_r))
      allocate(elec_rf_volr0_r(lsize_r))
      call mct_aVect_exportRattr(elec_rf_r_av, elec_rf_normalized_field, elec_rf_normalized_r)
      call mct_aVect_exportRattr(elec_rf_r_av, elec_rf_volr0_field, elec_rf_volr0_r)
  
      allocate(mfc_withd_normalized_r(lsize_r))
      allocate(mfc_withd_volr0_r(lsize_r))
      call mct_aVect_exportRattr(mfc_withd_r_av, mfc_withd_normalized_field, mfc_withd_normalized_r)
      call mct_aVect_exportRattr(mfc_withd_r_av, mfc_withd_volr0_field, mfc_withd_volr0_r)
  
      allocate(mfc_rf_normalized_r(lsize_r))
      allocate(mfc_rf_volr0_r(lsize_r))
      call mct_aVect_exportRattr(mfc_rf_r_av, mfc_rf_normalized_field, mfc_rf_normalized_r)
      call mct_aVect_exportRattr(mfc_rf_r_av, mfc_rf_volr0_field, mfc_rf_volr0_r)
  
      allocate(vmin_withd_normalized_r(lsize_r))
      allocate(vmin_withd_volr0_r(lsize_r))
      call mct_aVect_exportRattr(vmin_withd_r_av, vmin_withd_normalized_field, vmin_withd_normalized_r)
      call mct_aVect_exportRattr(vmin_withd_r_av, vmin_withd_volr0_field, vmin_withd_volr0_r)
  
      allocate(vmin_rf_normalized_r(lsize_r))
      allocate(vmin_rf_volr0_r(lsize_r))
      call mct_aVect_exportRattr(vmin_rf_r_av, vmin_rf_normalized_field, vmin_rf_normalized_r)
      call mct_aVect_exportRattr(vmin_rf_r_av, vmin_rf_volr0_field, vmin_rf_volr0_r)
  
      ! ------------------------------------------------------------------------
      ! Convert to a total sector water flux on the ROF grid, and put this in the l2r_rx
      ! attribute vector
      ! ------------------------------------------------------------------------
  
      allocate(dom_withd_flux_r(lsize_r))
      allocate(dom_rf_flux_r(lsize_r))
      allocate(liv_withd_flux_r(lsize_r))
      allocate(liv_rf_flux_r(lsize_r))
      allocate(elec_withd_flux_r(lsize_r))
      allocate(elec_rf_flux_r(lsize_r))
      allocate(mfc_withd_flux_r(lsize_r))
      allocate(mfc_rf_flux_r(lsize_r))
      allocate(min_withd_flux_r(lsize_r))
      allocate(min_rf_flux_r(lsize_r))
      do r = 1, lsize_r
         dom_withd_flux_r(r) = (dom_withd_normalized_r(r) * volr_r(r)) + dom_withd_volr0_r(r)
         dom_rf_flux_r(r) = (dom_withd_normalized_r(r) * volr_r(r)) + dom_withd_volr0_r(r)
  
         liv_withd_flux_r(r) = (liv_withd_normalized_r(r) * volr_r(r)) + liv_withd_volr0_r(r)
         liv_rf_flux_r(r) = (liv_withd_normalized_r(r) * volr_r(r)) + liv_withd_volr0_r(r)
  
         elec_withd_flux_r(r) = (elec_withd_normalized_r(r) * volr_r(r)) + elec_withd_volr0_r(r)
         elec_rf_flux_r(r) = (elec_withd_normalized_r(r) * volr_r(r)) + elec_withd_volr0_r(r)
  
         mfc_withd_flux_r(r) = (mfc_withd_normalized_r(r) * volr_r(r)) + mfc_withd_volr0_r(r)
         mfc_rf_flux_r(r) = (mfc_withd_normalized_r(r) * volr_r(r)) + mfc_withd_volr0_r(r)
  
         min_withd_flux_r(r) = (min_withd_normalized_r(r) * volr_r(r)) + min_withd_volr0_r(r)
         min_rf_flux_r(r) = (min_withd_normalized_r(r) * volr_r(r)) + min_withd_volr0_r(r)
      end do
  
      call mct_aVect_importRattr(l2r_r, dom_withd_flux_field, dom_withd_flux_r)
      call mct_aVect_importRattr(l2r_r, dom_rf_flux_field, dom_rf_flux_r)
      
      call mct_aVect_importRattr(l2r_r, liv_withd_flux_field, liv_withd_flux_r)
      call mct_aVect_importRattr(l2r_r, liv_rf_flux_field, liv_rf_flux_r)
  
      call mct_aVect_importRattr(l2r_r, elec_withd_flux_field, elec_withd_flux_r)
      call mct_aVect_importRattr(l2r_r, elec_rf_flux_field, elec_rf_flux_r)
  
      call mct_aVect_importRattr(l2r_r, mfc_withd_flux_field, mfc_withd_flux_r)
      call mct_aVect_importRattr(l2r_r, mfc_rf_flux_field, mfc_rf_flux_r)
  
      call mct_aVect_importRattr(l2r_r, min_withd_flux_field, min_withd_flux_r)
      call mct_aVect_importRattr(l2r_r, min_rf_flux_field, min_rf_flux_r)
  
  
  
      ! ------------------------------------------------------------------------
      ! Clean up
      ! ------------------------------------------------------------------------
  
      deallocate(volr_r)
      deallocate(volr_l)
  
      deallocate(dom_withd_flux_l)
      deallocate(dom_withd_flux_r)
      deallocate(dom_withd_normalized_l)
      deallocate(dom_withd_normalized_r)
      deallocate(dom_withd_volr0_l)
      deallocate(dom_withd_volr0_r)
      call mct_aVect_clean(dom_withd_l_av)
      call mct_aVect_clean(dom_withd_r_av)
      deallocate(dom_rf_flux_l)
      deallocate(dom_rf_flux_r)
      deallocate(dom_rf_normalized_l)
      deallocate(dom_rf_normalized_r)
      deallocate(dom_rf_volr0_l)
      deallocate(dom_rf_volr0_r)
      call mct_aVect_clean(dom_rf_l_av)
      call mct_aVect_clean(dom_rf_r_av)
  
      deallocate(elec_withd_flux_l)
      deallocate(elec_withd_flux_r)
      deallocate(elec_withd_normalized_l)
      deallocate(elec_withd_normalized_r)
      deallocate(elec_withd_volr0_l)
      deallocate(elec_withd_volr0_r)
      call mct_aVect_clean(elec_withd_l_av)
      call mct_aVect_clean(elec_withd_r_av)
      deallocate(elec_rf_flux_l)
      deallocate(elec_rf_flux_r)
      deallocate(elec_rf_normalized_l)
      deallocate(elec_rf_normalized_r)
      deallocate(elec_rf_volr0_l)
      deallocate(elec_rf_volr0_r)
      call mct_aVect_clean(elec_rf_l_av)
      call mct_aVect_clean(elec_rf_r_av)
  
      deallocate(liv_withd_flux_l)
      deallocate(liv_withd_flux_r)
      deallocate(liv_withd_normalized_l)
      deallocate(liv_withd_normalized_r)
      deallocate(liv_withd_volr0_l)
      deallocate(liv_withd_volr0_r)
      call mct_aVect_clean(liv_withd_l_av)
      call mct_aVect_clean(liv_withd_r_av)
      deallocate(liv_rf_flux_l)
      deallocate(liv_rf_flux_r)
      deallocate(liv_rf_normalized_l)
      deallocate(liv_rf_normalized_r)
      deallocate(liv_rf_volr0_l)
      deallocate(liv_rf_volr0_r)
      call mct_aVect_clean(liv_rf_l_av)
      call mct_aVect_clean(liv_rf_r_av)
  
      deallocate(mfc_withd_flux_l)
      deallocate(mfc_withd_flux_r)
      deallocate(mfc_withd_normalized_l)
      deallocate(mfc_withd_normalized_r)
      deallocate(mfc_withd_volr0_l)
      deallocate(mfc_withd_volr0_r)
      call mct_aVect_clean(mfc_withd_l_av)
      call mct_aVect_clean(mfc_withd_r_av)
      deallocate(mfc_rf_flux_l)
      deallocate(mfc_rf_flux_r)
      deallocate(mfc_rf_normalized_l)
      deallocate(mfc_rf_normalized_r)
      deallocate(mfc_rf_volr0_l)
      deallocate(mfc_rf_volr0_r)
      call mct_aVect_clean(mfc_rf_l_av)
      call mct_aVect_clean(mfc_rf_r_av)
  
      deallocate(min_withd_flux_l)
      deallocate(min_withd_flux_r)
      deallocate(min_withd_normalized_l)
      deallocate(min_withd_normalized_r)
      deallocate(min_withd_volr0_l)
      deallocate(min_withd_volr0_r)
      call mct_aVect_clean(min_withd_l_av)
      call mct_aVect_clean(min_withd_r_av)
      deallocate(min_rf_flux_l)
      deallocate(min_rf_flux_r)
      deallocate(min_rf_normalized_l)
      deallocate(min_rf_normalized_r)
      deallocate(min_rf_volr0_l)
      deallocate(min_rf_volr0_r)
      call mct_aVect_clean(min_rf_l_av)
      call mct_aVect_clean(min_rf_r_av)
  
    end subroutine map_lnd2rof_sectorwater
  
    subroutine map_rof2lnd_volr(volr_r, mapper_Fr2l, volr_l)
      !---------------------------------------------------------------
      ! Description
      ! Map volr from the rof grid to the lnd grid.
      !
      ! This is needed for the volr-normalization that is done in map_lnd2rof_sectorwater.
      !
      ! Note that this mapping is also done in the course of mapping all rof -> lnd fields
      ! in prep_lnd_calc_r2x_lx. However, we do this mapping ourselves here for two reasons:
      !
      ! (1) For the sake of this normalization, we change all volr < 0 values to 0; this is
      !     not done for the standard rof -> lnd mapping.
      !
      ! (2) It's possible that the driver sequencing would be changed such that this rof ->
      !     lnd mapping happens before the lnd -> rof mapping. If that happened, then volr_l
      !     (i.e., volr that has been mapped to the land grid by prep_lnd_calc_r2x_lx) would
      !     be inconsistent with volr_r, which would be a Bad Thing for the
      !     volr-normalizated mapping (this mapping would no longer be conservative). So we
      !     do the rof -> lnd remapping here to ensure we have a volr_l that is consistent
      !     with volr_r.
      !
      ! The pointer arguments to this routine should already be allocated to be the
      ! appropriate size.
      !
      ! Arguments
      real(r8), pointer, intent(in)    :: volr_r(:)   ! river volume on the rof grid (input)
      type(seq_map)    , intent(inout) :: mapper_Fr2l ! flux mapper for mapping rof -> lnd
      real(r8), pointer, intent(inout) :: volr_l(:)   ! river volume on the lnd grid (output) (technically intent(in) since intent gives the association status of a pointer, but given as intent(inout) to avoid confusion, since its data are modified)
      !
      ! Local variables
      integer :: lsize_r  ! number of rof points
      integer :: lsize_l  ! number of lnd points
      type(mct_avect) :: volr_r_av   ! temporary attribute vector holding volr on the rof grid
      type(mct_avect) :: volr_l_av   ! temporary attribute vector holding volr on the land grid
  
      ! This volr field name does not need to agree with the volr field name used in the
      ! 'real' attribute vectors
      character(len=*), parameter :: volr_field = 'volr'
      !---------------------------------------------------------------
  
      SHR_ASSERT_FL(associated(volr_r), sourcefile, __LINE__)
      SHR_ASSERT_FL(associated(volr_l), sourcefile, __LINE__)
  
      lsize_r = size(volr_r)
      lsize_l = size(volr_l)
  
      call mct_aVect_init(volr_r_av, rList = volr_field, lsize = lsize_r)
      call mct_aVect_importRattr(volr_r_av, volr_field, volr_r)
      call mct_aVect_init(volr_l_av, rList = volr_field, lsize = lsize_l)
  
      ! This mapping uses the same options as the standard rof -> lnd mapping done in
      ! prep_lnd_calc_r2x_lx. If that mapping ever changed (e.g., introducing an avwts_s
      ! argument), then it's *possible* that we'd want this mapping to change, too.
      call seq_map_map(mapper = mapper_Fr2l, &
           av_s = volr_r_av, &
           av_d = volr_l_av, &
           fldlist = volr_field, &
           norm = .true.)
  
      call mct_aVect_exportRattr(volr_l_av, volr_field, volr_l)
  
      call mct_aVect_clean(volr_r_av)
      call mct_aVect_clean(volr_l_av)
  
    end subroutine map_rof2lnd_volr
  
  end module map_lnd2rof_sectorwater_mod
  