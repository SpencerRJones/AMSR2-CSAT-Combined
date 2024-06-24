      MODULE define_csu1dvar

      implicit none
      save


     LOGICAL :: debug_mode = .TRUE.

      !--- these are shared variables, used by more than one 
      !--- subroutine and therefore commonly defined here
        
      character(12)  :: code_version = 'ver4.1'     
      character(12)  :: rtm_version  = 'MonoRTMv5.3' 
      character(100) :: salin_dir    = 'Aquarius_SSS/'   	! Aquarius salinity data directory
      character(100) :: bin_dir      = 'binary/'            ! Directory containing ERAI binary files     

      logical :: dontoutput      = .false. ! if true, code stops after OE step
      logical :: includeicecloud = .FALSE.  ! should match IWP being retrieved or not! (must also change nretvar below)

      !--- Variables to be modified for testing/running

      integer :: landlimit = 0 ! % land acceptable in FOV without skipping

      integer :: scstart =  1     ! starting scan number
      integer :: scend   =  50000 ! ending scan number
      integer :: pxstart =  1     ! starting pixel number
      integer :: pxend   =  5000  ! ending pixel number

      real    :: minlat =  -90.0 ! minimum latitude
      real    :: maxlat =   90.0 ! maximum latitude
      real    :: minlon = -180.0 ! minimum longitude
      real    :: maxlon =  180.0 ! maximum longitude

        INTEGER :: radar_stat = 0    !0 = Normal combined retrieval, 1 = radiometer-only retrieval, 2 = radar-only retrieval      

      !integer, parameter :: nz               = 16     ! number of layers in atmosphere -- cannot be modified!
      INTEGER, PARAMETER :: nz               = 30     ! modified: number of layers in the atmosphere
      INTEGER, PARAMETER :: nlevs            = 31     ! Number of height levels
      INTEGER            :: nlyrs
      INTEGER, PARAMETER :: nchans           = 10
      INTEGER, PARAMETER :: npas_freqs       = nchans
      INTEGER, PARAMETER :: nact_freqs       = 1
      INTEGER, PARAMETER :: nobs             = nchans + nz   !10 passive + 30 layers of reflectivity
      integer, parameter :: nbins            = 33     ! number of SST bins
      integer, parameter :: nang             = 30     ! number of EIA bins
      integer, parameter :: npc              =  3     ! number of water vapor EOFs/PCs to use (max=3)
      integer, parameter :: nn               =  20    ! max number of iterations to find convergence
      !integer, parameter :: nvar             =  7     ! number of parameters with potential to be retrieved
      !integer, parameter :: nretvar          =  4+npc ! number of parameters to retrieve (needs to be changed if includeicecloud is changed)
      !INTEGER, PARAMETER :: ncld             = nz   !Layers of cloud water
      !INTEGER, PARAMETER :: nplw             = nz   !Layers of rain water
      INTEGER, PARAMETER :: nvar             = 66 !Total number of retrieved variables
      INTEGER, PARAMETER :: nretvar          = nvar

      real :: conv_factor      = 20.0 ! value to divide nchan by to determine convergence (Rodgers Eq5.33)
      real :: chis1            =  1.0 ! chisq threshold for 'High Quality'
      real :: chisq_out_thresh =  100.0 ! if gt, oe_output set to missing

      ! Bounds on retrieved elements -- all that are capable of being retrieved!
      !                         EOF1  EOF2  EOF3  Wspd  log(LWP) log(IWP)    SST
      !                         ----  ----  ----  ----  -------- --------    ---
      !real :: x_min(nvar)! = (/ -4.0, -4.0, -4.0,  0.5,    -5.00,    -5.00, 270.0/)
      !real :: x_max(nvar)! = (/  4.0,  4.0,  4.0, 35.0,     0.12,      1.0, 307.0/) 

      integer,parameter :: knd1 = 1
      integer,parameter :: knd2 = 2
      integer,parameter :: knd8 = 8
      integer,parameter :: maxchans = 10
      integer,parameter :: nch = 10      !AMSR-2

      real               :: final_tpw              ! pass tpw between subroutines with this?
      character(len=256) :: orig_file              ! outfile for conversionmodule
      character(len=100) :: sy_file, sa_file
      character(len=10)  :: satellite, sensor
      character(len=20)  :: ppversion
      integer            :: nfreq 
      integer            :: nscans
      integer            :: npix
      integer            :: oepix, oelin
      integer            :: first_good_date        ! long integer string of YYYYMMDD'      
      integer(1)         :: chan_avail(maxchans)
      character(6)       :: freqs(maxchans)        ! whether ch exists or not

      !--- input variables

      real, ALLOCATABLE       :: plev_in(:,:)      ! pressure level profile
      real, allocatable       :: tbobs_in(:,:)    ! Observed Tb
      real, allocatable       :: eia_in(:,:)      ! earth incidence angle (pix,scan,chan)
      real, allocatable       :: alpha_in(:,:)	    ! scan angle
      real, allocatable       :: tlev_in(:,:)     ! temperature level profile
      !real, allocatable       :: hlev_in(:,:,:)     ! height level profile
      REAL, ALLOCATABLE       :: hlev_in(:)
      real, allocatable       :: rmix_in(:,:)     ! input pp wv mixing ratio profile
      REAL, ALLOCATABLE       :: refl_in(:,:)
      REAL, ALLOCATABLE       :: sst_in(:), rey_sst_in(:), rss_sst_in(:)
      REAL, ALLOCATABLE       :: t2m_in(:)
      REAL, ALLOCATABLE       :: d2m_in(:)         !2m dewpoint
      integer, allocatable    :: dtime_in(:,:)      ! date/time  (year, month, day, hour, minute, second)
      real, allocatable       :: lat_in(:)        ! pixel latitude
      real, allocatable       :: lon_in(:)        ! pixel longitude
      real, allocatable       :: sclat_in(:)        ! spacecraft latitude
      real, allocatable       :: sclon_in(:)        ! spacecraft longitude
      real, allocatable       :: scalt_in(:)        ! spacecraft altitude
      real, allocatable       :: scorient_in(:)     ! spacecraft orientation
      real, allocatable       :: slp_in(:)        ! sea level pressure
      real, allocatable       :: tpw_in(:)        ! input pp tpw
      real, allocatable       :: clwp_in(:), ciwp_in(:)
      real, allocatable       :: tsfc_in(:,:)       ! input pp surface skin temperature
      real, allocatable       :: wspd_in(:)       ! input pp surface wind speed
      integer*2, allocatable  :: wdir_in(:,:)       ! wind direction from analysis, deg from N
      integer*1, allocatable  :: land_in(:,:,:)     ! input pp land percent in FOV
      !integer*1, allocatable  :: sfc_in(:,:)        ! input surface code
      INTEGER, ALLOCATABLE    :: frzlvl_bin_in(:)
      REAL,    ALLOCATABLE    :: frzlvl_in(:)
      REAL,    ALLOCATABLE    :: cldbse_in(:)
      INTEGER, ALLOCATABLE    :: cldbse_bin_in(:)
      REAL,    ALLOCATABLE    :: lcl_in(:)
      INTEGER, ALLOCATABLE    :: qflg_in(:)
      INTEGER, ALLOCATABLE    :: sfc_in(:)
      INTEGER, ALLOCATABLE    :: modiscf_in(:,:)
      REAL, ALLOCATABLE       :: prcp_in(:)
      REAL, ALLOCATABLE       :: snrt_sfc_in(:)
      !--- output variables

      INTEGER, ALLOCATABLE :: cltop_bin(:), clbot_bin(:)

      real, allocatable       :: tbsim_out(:,:)   ! Simulated Tb
      real, allocatable       :: tbdif_out(:,:)   ! Observed minus simulated Tb difference
      real, allocatable       :: post_out(:,:)    ! posteriori error, stddev
      real, allocatable       :: saln_out(:)      ! salinity
      real, allocatable       :: tavg_out(:,:)    ! temperature profile output
      real, allocatable       :: rmix_out(:,:)    ! WV mixing ratio profile output
      real, allocatable       :: clwc_out(:,:)    ! cloud liquid water content
      real, allocatable       :: plwc_out(:,:)    ! precip liquid water content
      real, allocatable       :: sfcprcp_out(:)    ! precip rate 
      real, allocatable       :: rain_rate(:,:)
      real, allocatable       :: snow_rate(:,:)
      real, allocatable       :: tiwc_out(:,:)    ! total ice water content
      real, allocatable       :: incr_out(:,:)    ! increment info
      integer(1), allocatable :: iter_out(:)      ! number of iterations

      integer(1), allocatable :: lsmask(:,:)        ! landmask array
      integer(1), allocatable :: sfc_type(:)      ! surface type
      integer(1), allocatable :: lo_flag(:,:)       ! land/ocean flag
      !integer(1), allocatable :: sglint_in(:,:)     ! sun glint angle (degrees)
      real, allocatable       :: sst(:,:)           ! SST
      real, allocatable       :: sstgrid(:,:)       ! sst grid
      logical, allocatable    :: icegrid(:,:)       ! ice grid
      real,allocatable        :: freeqs(:)          ! channel frequencies
      real*8,allocatable      :: tai93time(:)       ! time in seconds since Jan 1 1993
      real,allocatable        :: pptcwv(:,:)        ! TPW from model, pass through

      real                    :: miss_flt = -9999.9  ! missing floating point value
      integer, parameter      :: miss_int = -9999    ! missing integer value
      integer, parameter      :: miss_byt = -99      ! missing byte value
      logical                 :: first = .TRUE.     ! indicates first time compute_tb is called

        integer, allocatable :: retr_flag_out(:)


      !--- OE variables passed back to retrieval

      real, allocatable    :: oe_output(:,:)
      real, allocatable    :: screen(:,:)
      real                 :: xr(nvar)          ! retrieval variables in OE
      real                 :: chisq             ! chi qquared
      real                 :: last_oeout(nvar)  ! used for a closer first guess 
      integer              :: loeop=-9,loeos=-9 ! save pix/scan for fg
      integer              :: run_count=0
      real                 :: avg_iter          ! averge number of iterations
        REAL, ALLOCATABLE :: z_a(:,:,:)
        REAL, ALLOCATABLE :: z_e(:,:,:)
        REAL, ALLOCATABLE :: atten(:)



      !--- variables for ERA-derived LUT table of means/sigmas

      character(len=90)      :: erafile
      integer*2, allocatable :: era_wm(:,:)
      integer*2, allocatable :: era_ws(:,:)
      integer*2              :: ssdex
      integer*2              :: agdex		! New - EIA index
      real, allocatable      :: mprof(:,:)
      real, allocatable      :: eofs(:,:,:)
      real, allocatable      :: peofs(:,:)
      real                   :: mrmp(nz)
      integer, parameter     :: gm=4            ! LUT grid multiplier -- set in LUT creation!
      integer, parameter     :: nlon=512*gm
      integer, parameter     :: nlat=256*gm     ! nmo=12
      real                   :: losize=0.703125 ! longitude bin size, need to match LUT creation!
      real                   :: lasize=0.701760 ! latitude bin size

      ! gsize is exact for lons, lats start at 89.463 and go by ~.701760

      real, allocatable :: sy(:,:) 
      real, allocatable    :: sy_dummy(:,:)
      real(8), allocatable :: sy_i(:,:)
      real, allocatable :: toffsets(:)
      real, allocatable :: s_toffsets(:,:)
      real, allocatable :: s_sy(:,:,:)

      REAL, ALLOCATABLE :: sa(:,:)

      ! variables for MonoRTM lookup table

      integer              :: nfreq_lut_pas, nfreq_lut_act
      integer              :: nchan_lut_pas, nchan_lut_act
      integer              :: npres_lut_pas, npres_lut_act
      integer              :: ntemp_lut_pas, ntemp_lut_act
      integer              :: nrmix_lut_pas, nrmix_lut_act

      integer, allocatable :: ifreq_lut_pas(:), ifreq_lut_act(:)
      integer, allocatable :: ipol_lut_pas(:)
      real, allocatable    :: freq_lut_pas(:), freq_lut_act(:)
      real, allocatable    :: pres_lut_pas(:), pres_lut_act(:)
      real, allocatable    :: temp_lut_pas(:), temp_lut_act(:)
      real, allocatable    :: rmix_lut_pas(:), rmix_lut_act(:)
      real, allocatable    :: kabs_lut(:,:,:,:)

      REAL, ALLOCATABLE :: kabs_lut_act(:,:,:,:)
      REAL, ALLOCATABLE :: kabs_lut_pas(:,:,:,:)


      integer              :: nfreq_lut
!      integer              :: nchan_lut
      integer              :: npres_lut
      integer              :: ntemp_lut
      integer              :: nrmix_lut

      integer, allocatable :: ifreq_lut(:)
      integer, allocatable :: ipol_lut(:)
      real, allocatable    :: freq_lut(:)
      real, allocatable    :: pres_lut(:)
      real, allocatable    :: temp_lut(:)
      real, allocatable    :: rmix_lut(:)
!      real, allocatable    :: kabs_lut(:,:,:,:)
      
      ! variables for cloud ice/snow lookup table
      
      integer              :: nfreq_icelut
      integer              :: ntemp_icelut
      integer              :: niwc_icelut
      integer              :: nrho_icelut

      integer, allocatable :: ifreq_icelut(:)
      real, allocatable    :: freq_icelut(:)
      real, allocatable    :: temp_icelut(:)
      real, allocatable    :: iwc_icelut(:)
      real, allocatable    :: rho_icelut(:)
      real, allocatable    :: ksca_icelut(:,:,:,:)
      real, allocatable    :: asca_icelut(:,:,:,:)
      real, allocatable    :: gsca_icelut(:,:,:,:)
      real, allocatable    :: pbck_icelut(:,:,:,:)

      real    :: stepsize_T
      real    :: stepsize_iwc
      real    :: stepsize_rho

      ! altitudes in [m], derived from mean geopotential at P levels

      real :: eof_lw_cov(3,nbins)
      real :: mr_sigmas(nbins,3)

      real, allocatable :: azim_in(:,:)   ! calculating satellite azimuthal angle

      real :: sss(720,360) ! sea surface salinity annual mean grid
        
      integer(4) :: granule
      logical    :: outdebug = .FALSE. ! if true will provide more output info for debugging

      !---DSD vars:

        !REAL, ALLOCATABLE :: N0_lut(:)
        REAL              :: N0_lut
        REAL, ALLOCATABLE :: mu_lut(:)
        REAL, ALLOCATABLE :: lamb_lut(:)
        !REAL, ALLOCATABLE :: gamma_lwc_lut(:,:,:)
        REAL, ALLOCATABLE :: gamma_lwc_lut(:,:)






      END MODULE define_csu1dvar
