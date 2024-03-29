module shr_drydep_mod

  !========================================================================
  ! Module for handling dry depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !
  ! !REVISION HISTORY:
  !     2008-Nov-12 - F. Vitt - creation.
  !     2009-Feb-19 - E. Kluzek - merge shr_drydep_tables module in.
  !     2009-Feb-20 - E. Kluzek - use shr_ coding standards, and check for namelist file.
  !     2009-Feb-20 - E. Kluzek - Put _r8 on all constants, remove namelist read out.
  !     2009-Mar-23 - F. Vitt - Some corrections/cleanup and addition of drydep_method.
  !     2009-Mar-27 - E. Kluzek - Get description and units from J.F. Lamarque.
  !========================================================================

  ! !USES:

  use shr_sys_mod,   only : shr_sys_abort
  use shr_log_mod,   only : s_loglev  => shr_log_Level
  use shr_kind_mod,  only : r8 => shr_kind_r8, CS => SHR_KIND_CS
  use shr_const_mod, only : SHR_CONST_MWWV

  implicit none
  save

  private

  ! !PUBLIC MEMBER FUNCTIONS

  public :: shr_drydep_readnl       ! Read namelist
  public :: shr_drydep_init         ! Initialization of drydep data
  public :: shr_drydep_setHCoeff    ! Calculate Henry's law coefficients

  ! !PRIVATE ARRAY SIZES

  integer, private, parameter :: maxspc = 210              ! Maximum number of species
  integer, public,  protected :: n_species_table           ! Number of species to work with
  integer, private, parameter :: NSeas = 5                 ! Number of seasons
  integer, public,  parameter :: NLUse = 11                ! Number of land-use types
  integer, private, protected :: NHen

  ! !PUBLIC DATA MEMBERS:

  real(r8), public, parameter :: ph     = 1.e-5_r8         ! measure of the acidity (dimensionless)

  integer, public, protected  :: n_drydep = 0                         ! Number in drypdep list
  character(len=32), public, protected :: drydep_list(maxspc) = ''   ! List of dry-dep species

  character(len=CS), public, protected :: drydep_fields_token = ''   ! First drydep fields token

  real(r8), public, allocatable, protected :: foxd(:)      ! reactivity factor for oxidation (dimensioness)
  real(r8), public, allocatable, protected :: drat(:)      ! ratio of molecular diffusivity (D_H2O/D_species; dimensionless)
  integer,  public, allocatable, protected :: mapping(:)   ! mapping to species table
  ! --- Indices for each species ---
  integer,  public, protected :: h2_ndx, ch4_ndx, co_ndx, pan_ndx, mpan_ndx, so2_ndx, o3_ndx, o3a_ndx, xpan_ndx

  !---------------------------------------------------------------------------
  ! Table 1 from Wesely, Atmos. Environment, 1989, p1293
  ! Table 2 from Sheih, microfiche PB86-218104 and Walcek, Atmos.  Environment, 1986, p949
  ! Table 3-5 compiled by P. Hess
  !
  ! index #1 : season
  !           1 -> midsummer with lush vegetation
  !           2 -> autumn with unharvested cropland
  !           3 -> late autumn after frost, no snow
  !           4 -> winter, snow on ground, and subfreezing
  !           5 -> transitional spring with partially green short annuals
  !
  ! index #2 : landuse type
  !           1 -> urban land
  !           2 -> agricultural land
  !           3 -> range land
  !           4 -> deciduous forest
  !           5 -> coniferous forest
  !           6 -> mixed forest including wetland
  !           7 -> water, both salt and fresh
  !           8 -> barren land, mostly desert
  !           9 -> nonforested wetland
  !           10 -> mixed agricultural and range land
  !           11 -> rocky open areas with low growing shrubs
  !
  ! JFL August 2000
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! table to parameterize the impact of soil moisture on the deposition of H2 and
  ! CO on soils (from Sanderson et al., J. Atmos. Chem., 46, 15-28, 2003).
  !---------------------------------------------------------------------------

  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_a(NLUse) = &
       (/  0.000_r8,  0.000_r8, 0.270_r8,  0.000_r8,  0.000_r8,  &
       0.000_r8,  0.000_r8, 0.000_r8,  0.000_r8,  0.000_r8, 0.000_r8/)
  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_b(NLUse) = &
       (/  0.000_r8,-41.390_r8, -0.472_r8,-41.900_r8,-41.900_r8,  &
       -41.900_r8,  0.000_r8,  0.000_r8,  0.000_r8,-41.390_r8,  0.000_r8/)
  !--- deposition of h2 and CO on soils ---
  real(r8), parameter, public :: h2_c(NLUse) = &
       (/  0.000_r8, 16.850_r8, 1.235_r8, 19.700_r8, 19.700_r8, &
       19.700_r8,  0.000_r8, 0.000_r8,  0.000_r8, 17.700_r8, 1.000_r8/)

  !--- deposition of h2 and CO on soils
  !
  !--- ri:   Richardson number                      (dimensionless)
  !--- rlu:  Resistance of leaves in upper canopy   (s.m-1)
  !--- rac:  Aerodynamic resistance to lower canopy (s.m-1)
  !--- rgss: Ground surface resistance for SO2      (s.m-1)
  !--- rgso: Ground surface resistance for O3       (s.m-1)
  !--- rcls: Lower canopy resistance for SO2        (s.m-1)
  !--- rclo: Lower canopy resistance for O3         (s.m-1)
  !
  real(r8), public, protected, dimension(NSeas,NLUse) :: ri, rlu, rac, rgss, rgso, rcls, rclo

  data ri  (1,1:NLUse) &
       /1.e36_r8,  60._r8, 120._r8,  70._r8, 130._r8, 100._r8,1.e36_r8,1.e36_r8,  80._r8, 100._r8, 150._r8/
  data rlu (1,1:NLUse) &
       /1.e36_r8,2000._r8,2000._r8,2000._r8,2000._r8,2000._r8,1.e36_r8,1.e36_r8,2500._r8,2000._r8,4000._r8/
  data rac (1,1:NLUse) &
       / 100._r8, 200._r8, 100._r8,2000._r8,2000._r8,2000._r8,   0._r8,   0._r8, 300._r8, 150._r8, 200._r8/
  data rgss(1,1:NLUse) &
       / 400._r8, 150._r8, 350._r8, 500._r8, 500._r8, 100._r8,   0._r8,1000._r8,  0._r8, 220._r8, 400._r8/
  data rgso(1,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(1,1:NLUse) &
       /1.e36_r8,2000._r8,2000._r8,2000._r8,2000._r8,2000._r8,1.e36_r8,1.e36_r8,2500._r8,2000._r8,4000._r8/
  data rclo(1,1:NLUse) &
       /1.e36_r8,1000._r8,1000._r8,1000._r8,1000._r8,1000._r8,1.e36_r8,1.e36_r8,1000._r8,1000._r8,1000._r8/

  data ri  (2,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 250._r8, 500._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (2,1:NLUse) &
       /1.e36_r8,9000._r8,9000._r8,9000._r8,4000._r8,8000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (2,1:NLUse) &
       / 100._r8, 150._r8, 100._r8,1500._r8,2000._r8,1700._r8,   0._r8,   0._r8, 200._r8, 120._r8, 140._r8/
  data rgss(2,1:NLUse) &
       / 400._r8, 200._r8, 350._r8, 500._r8, 500._r8, 100._r8,   0._r8,1000._r8,   0._r8, 300._r8, 400._r8/
  data rgso(2,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8, 800._r8, 180._r8, 200._r8/
  data rcls(2,1:NLUse) &
       /1.e36_r8,9000._r8,9000._r8,9000._r8,2000._r8,4000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rclo(2,1:NLUse) &
       /1.e36_r8, 400._r8, 400._r8, 400._r8,1000._r8, 600._r8,1.e36_r8,1.e36_r8, 400._r8, 400._r8, 400._r8/

  data ri  (3,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 250._r8, 500._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (3,1:NLUse) &
       /1.e36_r8,1.e36_r8,9000._r8,9000._r8,4000._r8,8000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (3,1:NLUse) &
       / 100._r8,  10._r8, 100._r8,1000._r8,2000._r8,1500._r8,   0._r8,   0._r8, 100._r8, 50._r8, 120._r8/
  data rgss(3,1:NLUse) &
       / 400._r8, 150._r8, 350._r8, 500._r8, 500._r8, 200._r8,   0._r8,1000._r8,   0._r8, 200._r8, 400._r8/
  data rgso(3,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(3,1:NLUse) &
       /1.e36_r8,1.e36_r8,9000._r8,9000._r8,3000._r8,6000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rclo(3,1:NLUse) &
       /1.e36_r8,1000._r8, 400._r8, 400._r8,1000._r8, 600._r8,1.e36_r8,1.e36_r8, 800._r8, 600._r8, 600._r8/

  data ri  (4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8, 400._r8, 800._r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8/
  data rlu (4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,1.e36_r8,6000._r8,9000._r8,1.e36_r8,1.e36_r8,9000._r8,9000._r8,9000._r8/
  data rac (4,1:NLUse) &
       / 100._r8,  10._r8,  10._r8,1000._r8,2000._r8,1500._r8,   0._r8,   0._r8,  50._r8,  10._r8,  50._r8/
  data rgss(4,1:NLUse) &
       / 100._r8, 100._r8, 100._r8, 100._r8, 100._r8, 100._r8,   0._r8,1000._r8, 100._r8, 100._r8,  50._r8/
  data rgso(4,1:NLUse) &
       / 600._r8,3500._r8,3500._r8,3500._r8,3500._r8,3500._r8,2000._r8, 400._r8,3500._r8,3500._r8,3500._r8/
  data rcls(4,1:NLUse) &
       /1.e36_r8,1.e36_r8,1.e36_r8,9000._r8, 200._r8, 400._r8,1.e36_r8,1.e36_r8,9000._r8,1.e36_r8,9000._r8/
  data rclo(4,1:NLUse) &
       /1.e36_r8,1000._r8,1000._r8, 400._r8,1500._r8, 600._r8,1.e36_r8,1.e36_r8, 800._r8,1000._r8, 800._r8/

  data ri  (5,1:NLUse) &
       /1.e36_r8, 120._r8, 240._r8, 140._r8, 250._r8, 190._r8,1.e36_r8,1.e36_r8, 160._r8, 200._r8, 300._r8/
  data rlu (5,1:NLUse) &
       /1.e36_r8,4000._r8,4000._r8,4000._r8,2000._r8,3000._r8,1.e36_r8,1.e36_r8,4000._r8,4000._r8,8000._r8/
  data rac (5,1:NLUse) &
       / 100._r8,  50._r8,  80._r8,1200._r8,2000._r8,1500._r8,   0._r8,   0._r8, 200._r8, 60._r8, 120._r8/
  data rgss(5,1:NLUse) &
       / 500._r8, 150._r8, 350._r8, 500._r8, 500._r8, 200._r8,   0._r8,1000._r8,   0._r8, 250._r8, 400._r8/
  data rgso(5,1:NLUse) &
       / 300._r8, 150._r8, 200._r8, 200._r8, 200._r8, 300._r8,2000._r8, 400._r8,1000._r8, 180._r8, 200._r8/
  data rcls(5,1:NLUse) &
       /1.e36_r8,4000._r8,4000._r8,4000._r8,2000._r8,3000._r8,1.e36_r8,1.e36_r8,4000._r8,4000._r8,8000._r8/
  data rclo(5,1:NLUse) &
       /1.e36_r8,1000._r8, 500._r8, 500._r8,1500._r8, 700._r8,1.e36_r8,1.e36_r8, 600._r8, 800._r8, 800._r8/

  !---------------------------------------------------------------------------
  !         ... roughness length
  !---------------------------------------------------------------------------
  real(r8), public, protected, dimension(NSeas,NLUse) :: z0

  data z0  (1,1:NLUse) &
       /1.000_r8,0.250_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.150_r8,0.100_r8,0.100_r8/
  data z0  (2,1:NLUse) &
       /1.000_r8,0.100_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.100_r8,0.080_r8,0.080_r8/
  data z0  (3,1:NLUse) &
       /1.000_r8,0.005_r8,0.050_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.100_r8,0.020_r8,0.060_r8/
  data z0  (4,1:NLUse) &
       /1.000_r8,0.001_r8,0.001_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.001_r8,0.001_r8,0.040_r8/
  data z0  (5,1:NLUse) &
       /1.000_r8,0.030_r8,0.020_r8,1.000_r8,1.000_r8,1.000_r8,0.0006_r8,0.002_r8,0.010_r8,0.030_r8,0.060_r8/

  !---------------------------------------------------------------------------
  ! public chemical data
  !---------------------------------------------------------------------------

  !--- data for foxd (reactivity factor for oxidation) ----
  real(r8), public, protected, allocatable :: dfoxd(:)

  ! PRIVATE DATA:

  Interface shr_drydep_setHCoeff
     Module Procedure set_hcoeff_scalar
     Module Procedure set_hcoeff_vector
  End Interface shr_drydep_setHCoeff

  real(r8), private, parameter :: small_value = 1.e-36_r8          !--- smallest value to use ---

  !---------------------------------------------------------------------------
  ! private chemical data
  !---------------------------------------------------------------------------

  !--- Names of species that can work with ---
  character(len=16), public, protected, allocatable :: species_name_table(:)

  !--- data for effective Henry's Law coefficient ---
  real(r8), public, protected, allocatable :: dheff(:,:)

  real(r8), parameter :: wh2o = SHR_CONST_MWWV
  real(r8), allocatable :: mol_wgts(:)

  character(len=500) :: dep_data_file = 'NONE' ! complete file path

  !===============================================================================
CONTAINS
  !===============================================================================

  !====================================================================================

  subroutine shr_drydep_readnl(NLFilename, ID, shr_drydep_fields)

    !========================================================================
    ! reads drydep_inparm namelist and sets up CCSM driver list of fields for
    ! land-atmosphere communications.
    !
    ! !REVISION HISTORY:
    !  2009-Feb-20 - E. Kluzek - Separate out as subroutine from previous input_init
    !========================================================================

    use shr_file_mod,only : shr_file_getUnit, shr_file_freeUnit
    use shr_log_mod, only : s_logunit => shr_log_Unit
    use seq_comm_mct,only : seq_comm_iamroot, seq_comm_setptrs
    use shr_mpi_mod, only : shr_mpi_bcast
    use shr_nl_mod, only : shr_nl_find_group_name
    implicit none

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer         , intent(in)  :: ID         ! seq_comm ID
    character(len=*), intent(out) :: shr_drydep_fields

    !----- local -----
    integer :: i                ! Indices
    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    character(len=8) :: token   ! dry dep field name to add
    integer :: mpicom           ! MPI communicator

    !----- formats -----
    character(*),parameter :: subName = '(shr_drydep_read) '
    character(*),parameter :: F00   = "('(shr_drydep_read) ',8a)"
    character(*),parameter :: FI1   = "('(shr_drydep_init) ',a,I2)"

    namelist /drydep_inparm/ drydep_list, dep_data_file

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the drydep field list to pass
    ! First check if file exists and if not, n_drydep will be zero
    !-----------------------------------------------------------------------------

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0  )then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if
    call seq_comm_setptrs(ID,mpicom=mpicom)
    if (seq_comm_iamroot(ID)) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          unitn = shr_file_getUnit()
          open( unitn, file=trim(NLFilename), status='old' )
          if ( s_loglev > 0 ) write(s_logunit,F00) &
               'Read in drydep_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'drydep_inparm', ierr)
          if (ierr == 0) then
             ierr = 1
             do while ( ierr /= 0 )
                read(unitn, drydep_inparm, iostat=ierr)
                if (ierr < 0) then
                   call shr_sys_abort( subName//'ERROR: encountered end-of-file on namelist read' )
                endif
             end do
          else
             write(s_logunit,*) 'shr_drydep_read:  no drydep_inparm namelist found in ',NLFilename
          endif
          close( unitn )
          call shr_file_freeUnit( unitn )
       end if
    end if
    call shr_mpi_bcast( drydep_list, mpicom )
    call shr_mpi_bcast( dep_data_file, mpicom )

    n_drydep = 0

    !--- Loop over species to fill list of fields to communicate for drydep ---
    shr_drydep_fields = ' '
    do i=1,maxspc
       if ( len_trim(drydep_list(i))==0 ) exit
       write(token,333) i
       shr_drydep_fields = trim(shr_drydep_fields)//':'//trim(token)
       if ( i == 1 ) then
          shr_drydep_fields = trim(token)
          drydep_fields_token = trim(token)
       endif
       n_drydep = n_drydep+1
    enddo

    if ( s_loglev > 0 ) then
       if ( n_drydep == 0 )then
          write(s_logunit,F00) 'No dry deposition fields will be transfered'
       else
          write(s_logunit,FI1) 'Number of dry deposition fields transfered is ', &
               n_drydep
       end if
    end if

    ! Need to explicitly add Sl_ based on naming convention
333 format ('Sl_dd',i3.3)

  end subroutine shr_drydep_readnl

  !====================================================================================

  subroutine shr_drydep_init( )

    !========================================================================
    ! Initialization of dry deposition fields
    ! reads drydep_inparm namelist and sets up CCSM driver list of fields for
    ! land-atmosphere communications.
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !  2009-Feb-20 - E. Kluzek - Check for existance of file if not return, set n_drydep=0
    !  2009-Feb-20 - E. Kluzek - Move namelist read to separate subroutine
    !========================================================================

    use shr_log_mod, only : s_logunit => shr_log_Unit
    use shr_infnan_mod, only: shr_infnan_posinf, assignment(=)

    use shr_pio_mod,  only: shr_pio_getiosys, shr_pio_getiotype
    use seq_comm_mct
    use pio

    implicit none

    !----- local -----
    integer :: i, l                      ! Indices
    character(len=32) :: test_name       ! field test name
    !----- formats -----
    character(*),parameter :: subName = '(shr_drydep_init) '
    character(*),parameter :: F00   = "('(shr_drydep_init) ',8a)"

    type(file_desc_t) :: fileid

    character(len=seq_comm_namelen) :: cpl_name
    type(iosystem_desc_t) , pointer :: cpl_io_subsystem
    integer                         :: cpl_pio_iotype
    integer :: rcode, dimid, varid

    if (dep_data_file=='NONE' .or. len_trim(dep_data_file)==0) return

    ! read data from input file
    cpl_name         =  seq_comm_name(CPLID)
    cpl_io_subsystem => shr_pio_getiosys(cpl_name)
    cpl_pio_iotype   =  shr_pio_getiotype(cpl_name)

    rcode = pio_openfile(cpl_io_subsystem, fileid, cpl_pio_iotype, trim(dep_data_file),pio_nowrite)

    rcode = pio_inq_dimid(fileid,'n_species_table',dimid)
    rcode = pio_inquire_dimension(fileid,dimid,len=n_species_table)
    rcode = pio_inq_dimid(fileid,'NHen',dimid)
    rcode = pio_inquire_dimension(fileid,dimid,len=nHen)

    allocate( mol_wgts(n_species_table) )
    allocate( dfoxd(n_species_table) )
    allocate( species_name_table(n_species_table) )
    allocate( dheff(nhen,n_species_table))

    rcode = pio_inq_varid(fileid,'mol_wghts',varid)
    rcode = pio_get_var(fileid,varid,mol_wgts)

    rcode = pio_inq_varid(fileid,'dfoxd',varid)
    rcode = pio_get_var(fileid,varid,dfoxd)

    rcode = pio_inq_varid(fileid,'species_name_table',varid)
    rcode = pio_get_var(fileid,varid,species_name_table)

    rcode = pio_inq_varid(fileid,'dheff',varid)
    rcode = pio_get_var(fileid,varid,dheff)

    call pio_closefile(fileid)

    !-----------------------------------------------------------------------------
    ! Allocate and fill foxd, drat and mapping as well as species indices
    !-----------------------------------------------------------------------------

    if ( n_drydep > 0 ) then

       allocate( foxd(n_drydep) )
       allocate( drat(n_drydep) )
       allocate( mapping(n_drydep) )

       ! This initializes these variables to infinity.
       foxd = shr_infnan_posinf
       drat = shr_infnan_posinf

       mapping(:) = 0

    end if

    h2_ndx=-1; ch4_ndx=-1; co_ndx=-1; mpan_ndx = -1; pan_ndx = -1; so2_ndx=-1; o3_ndx=-1; xpan_ndx=-1

    !--- Loop over drydep species that need to be worked with ---
    do i=1,n_drydep
       if ( len_trim(drydep_list(i))==0 ) exit

       test_name = drydep_list(i)

       if( trim(test_name) == 'O3' ) then
          test_name = 'OX'
       end if

       !--- Figure out if species maps to a species in the species table ---
       do l = 1,n_species_table
          if(  trim( test_name ) == trim( species_name_table(l) ) ) then
             mapping(i)  = l
             exit
          end if
       end do

       !--- If it doesn't map to a species in the species table find species close enough ---
       if( mapping(i) < 1 ) then
          select case( trim(test_name) )
          case( 'O3S', 'O3INERT' )
             test_name = 'OX'
          case( 'Pb' )
             test_name = 'HNO3'
          case( 'SOGM','SOGI','SOGT','SOGB','SOGX' )
             test_name = 'CH3OOH'
          case( 'SOA', 'SO4', 'CB1', 'CB2', 'OC1', 'OC2', 'NH4', 'SA1', 'SA2', 'SA3', 'SA4' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'SOAM', 'SOAI', 'SOAT', 'SOAB', 'SOAX' )
             test_name = 'OX'  ! this is just a place holder. values are explicitly set below
          case( 'SOAGbb0' )
             test_name = 'SOAGff0'
          case( 'SOAGbb1' )
             test_name = 'SOAGff1'
          case( 'SOAGbb2' )
             test_name = 'SOAGff2'
          case( 'SOAGbb3' )
             test_name = 'SOAGff3'
          case( 'SOAGbb4' )
             test_name = 'SOAGff4'
          case( 'O3A' )
             test_name = 'OX'
          case( 'XMPAN' )
             test_name = 'MPAN'
          case( 'XPAN' )
             test_name = 'PAN'
          case( 'XNO' )
             test_name = 'NO'
          case( 'XNO2' )
             test_name = 'NO2'
          case( 'XHNO3' )
             test_name = 'HNO3'
          case( 'XONIT' )
             test_name = 'ONIT'
          case( 'XONITR' )
             test_name = 'ONITR'
          case( 'XHO2NO2')
             test_name = 'HO2NO2'
          case( 'XNH4NO3' )
             test_name = 'HNO3'
          case( 'NH4NO3' )
             test_name = 'HNO3'
          case default
             test_name = 'blank'
          end select

          !--- If found a match check the species table again ---
          if( trim(test_name) /= 'blank' ) then
             do l = 1,n_species_table
                if( trim( test_name ) == trim( species_name_table(l) ) ) then
                   mapping(i)  = l
                   exit
                end if
             end do
          else
             if ( s_loglev > 0 ) write(s_logunit,F00) trim(drydep_list(i)), &
                  ' not in tables; will have dep vel = 0'
             call shr_sys_abort( subName//': '//trim(drydep_list(i))//' is not in tables' )
          end if
       end if

       !--- Figure out the specific species indices ---
       if ( trim(drydep_list(i)) == 'H2' )   h2_ndx   = i
       if ( trim(drydep_list(i)) == 'CO' )   co_ndx   = i
       if ( trim(drydep_list(i)) == 'CH4' )  ch4_ndx  = i
       if ( trim(drydep_list(i)) == 'MPAN' ) mpan_ndx = i
       if ( trim(drydep_list(i)) == 'PAN' )  pan_ndx  = i
       if ( trim(drydep_list(i)) == 'SO2' )  so2_ndx  = i
       if ( trim(drydep_list(i)) == 'OX' .or. trim(drydep_list(i)) == 'O3' ) o3_ndx  = i
       if ( trim(drydep_list(i)) == 'O3A' ) o3a_ndx  = i
       if ( trim(drydep_list(i)) == 'XPAN' ) xpan_ndx = i

       if( mapping(i) > 0) then
          l = mapping(i)
          foxd(i) = dfoxd(l)
          drat(i) = sqrt(mol_wgts(l)/wh2o)
       endif

    enddo

    where( rgss < 1._r8 )
       rgss = 1._r8
    endwhere

    where( rac < small_value)
       rac = small_value
    endwhere

  end subroutine shr_drydep_init

  !====================================================================================

  subroutine set_hcoeff_scalar( sfc_temp, heff )

    !========================================================================
    ! Interface to shr_drydep_setHCoeff when input is scalar
    ! wrapper routine used when surface temperature is a scalar (single column) rather
    ! than an array (multiple columns).
    !
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !========================================================================

    implicit none

    real(r8), intent(in)     :: sfc_temp         ! Input surface temperature
    real(r8), intent(out)    :: heff(n_drydep)   ! Output Henry's law coefficients

    !----- local -----
    real(r8) :: sfc_temp_tmp(1)    ! surface temp

    sfc_temp_tmp(:) = sfc_temp
    call set_hcoeff_vector( 1, sfc_temp_tmp, heff(:n_drydep) )

  end subroutine set_hcoeff_scalar

  !====================================================================================

  subroutine set_hcoeff_vector( ncol, sfc_temp, heff )

    !========================================================================
    ! Interface to shr_drydep_setHCoeff when input is vector
    ! sets dry depositions coefficients -- used by both land and atmosphere models
    ! !REVISION HISTORY:
    !  2008-Nov-12 - F. Vitt - first version
    !========================================================================

    use shr_log_mod, only : s_logunit => shr_log_Unit

    implicit none

    integer, intent(in)      :: ncol                  ! Input size of surface-temp vector
    real(r8), intent(in)     :: sfc_temp(ncol)        ! Surface temperature
    real(r8), intent(out)    :: heff(ncol,n_drydep)   ! Henry's law coefficients

    !----- local -----
    real(r8), parameter :: t0     = 298._r8    ! Standard Temperature
    real(r8), parameter :: ph_inv = 1._r8/ph   ! Inverse of PH
    integer  :: m, l           ! indices
    real(r8) :: e298           ! Henry's law coefficient @ standard temperature (298K)
    real(r8) :: dhr            ! temperature dependence of Henry's law coefficient
    real(r8) :: dk1s(ncol)     ! DK Work array 1
    real(r8) :: dk2s(ncol)     ! DK Work array 2
    real(r8) :: wrk(ncol)      ! Work array

    !----- formats -----
    character(*),parameter :: subName = '(shr_drydep_set_hcoeff) '
    character(*),parameter :: F00   = "('(shr_drydep_set_hcoeff) ',8a)"

    !-------------------------------------------------------------------------------
    ! notes:
    !-------------------------------------------------------------------------------

    wrk(:) = (t0 - sfc_temp(:))/(t0*sfc_temp(:))
    do m = 1,n_drydep
       l    = mapping(m)
       e298 = dheff(1,l)
       dhr  = dheff(2,l)
       heff(:,m) = e298*exp( dhr*wrk(:) )
       !--- Calculate coefficients based on the drydep tables ---
       if( dheff(3,l) /= 0._r8 .and. dheff(5,l) == 0._r8 ) then
          e298 = dheff(3,l)
          dhr  = dheff(4,l)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          where( heff(:,m) /= 0._r8 )
             heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph_inv)
          elsewhere
             heff(:,m) = dk1s(:)*ph_inv
          endwhere
       end if
       !--- For coefficients that are non-zero AND CO2 or NH3 handle things this way ---
       if( dheff(5,l) /= 0._r8 ) then
          if( trim( drydep_list(m) ) == 'CO2' .or. trim( drydep_list(m) ) == 'NH3' .or. trim( drydep_list(m) ) == 'SO2' ) then
             e298 = dheff(3,l)
             dhr  = dheff(4,l)
             dk1s(:) = e298*exp( dhr*wrk(:) )
             e298 = dheff(5,l)
             dhr  = dheff(6,l)
             dk2s(:) = e298*exp( dhr*wrk(:) )
             !--- For Carbon dioxide ---
             if( trim(drydep_list(m)) == 'CO2'.or. trim( drydep_list(m) ) == 'SO2' ) then
                heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph_inv*(1._r8 + dk2s(:)*ph_inv))
                !--- For NH3 ---
             else if( trim( drydep_list(m) ) == 'NH3' ) then
                heff(:,m) = heff(:,m)*(1._r8 + dk1s(:)*ph/dk2s(:))
                !--- This can't happen ---
             else
                write(s_logunit,F00) 'Bad species ',drydep_list(m)
                call shr_sys_abort( subName//'ERROR: in assigning coefficients' )
             end if
          end if
       end if
    end do

  end subroutine set_hcoeff_vector

  !===============================================================================

end module shr_drydep_mod
