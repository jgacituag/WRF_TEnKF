PROGRAM make_synth_obs
!=======================================================================
! Synthetic radar reflectivity on model grid points, for WRF OSSEs.
!
! WRF_TEnKF addition (J. Gacitua). It is built on the WRF-LETKF modules of
! T. Miyoshi and J. Ruiz, so that the truth and the ensemble see exactly
! the observation operator of the analysis (Trans_XtoY -> calc_ref_vr).
!
! 1. Every rank reads ensemble members (ens_prefix//'NNNNN', NNNNN = 1..nens)
!    and computes the reflectivity at every grid point; the number of members
!    with echo (Z > 10**(dbz_min/10)) is summed over ranks: n_active.
! 2. Rank 0 reads the truth and computes h(x_true) in dBZ, censored at
!    dbz_min. It draws one N(0,1) number per candidate point, in a fixed
!    order (SFMT, fixed seed), so the noise does not depend on the ensemble:
!      y = max( h(x_true) + obs_err * eps , dbz_min )
! 3. Rank 0 writes the candidates with n_active >= min_active as a radar file
!    meant for radar_obs_lonlatz = .true.:
!      header : lon, lat, z of a nominal radar (domain centre, 0 m)
!      records: id_reflectivity_obs, lon, lat, z [m asl], Z [mm6 m-3],
!               obs_err [dBZ], radar_wavelength [cm]
!    with Z = 10**(y/10) if y > dbz_min and Z = 0 otherwise (the LETKF maps
!    Z <= 10**(dbz_min/10) to dbz_min). truth_obs_file holds the same points
!    with the noise-free h(x_true), for verification.
!
! Candidates: i = 3..nlon-2, j = 3..nlat-2, k = kmin..nlev-2 with z <= zmax,
! thinned by skip_ij and skip_k. z = PH/g of the truth, the height the
! LETKF uses to place radar observations (z2k_fast). Positions come from
! ij_to_latlon with the model projection, the inverse of what the LETKF
! applies, and the margins keep every point inside the range it accepts.
!
! Input : synth_obs.namelist, group &synth_obs (defaults below), the truth
!         file and the ensemble files.
! Output: obs_file, truth_obs_file, and a summary on unit 6 of rank 0.
! Run with an unlimited stack (ulimit -s unlimited): Trans_XtoY declares
! three automatic arrays of the size of the model grid.
!=======================================================================
  USE common
  USE common_mpi
  USE common_wrf
  USE common_namelist
  USE common_obs_wrf
  USE map_utils
  IMPLICIT NONE
  REAL(r_dble), EXTERNAL :: genrand_res53   !SFMT.f90 (external procedures)
  EXTERNAL :: init_gen_rand

  !Namelist
  CHARACTER(256) :: truth_file     = 'truth'
  CHARACTER(256) :: ens_prefix     = 'gs01'
  INTEGER        :: nens           = 59
  REAL(r_size)   :: obs_err        = 2.2360679775d0  ![dBZ], R = 5 dBZ**2
  REAL(r_size)   :: dbz_min        = 0.0d0           ![dBZ], detection limit
  INTEGER        :: seed           = 60
  REAL(r_size)   :: zmax           = 16000.0d0       ![m asl]
  INTEGER        :: kmin           = 2
  INTEGER        :: skip_ij        = 1
  INTEGER        :: skip_k         = 1
  INTEGER        :: min_active     = 1
  REAL(r_size)   :: radar_wavelength = 10.0d0        ![cm], 10 -> S band (Tong and Xue 2006)
  CHARACTER(256) :: obs_file       = 'rad0101.dat'
  CHARACTER(256) :: truth_obs_file = 'truth_rad0101.dat'
  NAMELIST / synth_obs / truth_file, ens_prefix, nens, obs_err, dbz_min, seed, zmax, kmin, &
                         skip_ij, skip_k, min_active, radar_wavelength, obs_file, truth_obs_file

  REAL(r_size), ALLOCATABLE :: v3d(:,:,:,:), v2d(:,:,:), zt(:,:,:), ht(:,:,:)
  INTEGER(4),   ALLOCATABLE :: nact_loc(:,:,:), nact(:,:,:)
  REAL(r_size) :: elm, typ, odat, zero, href, eps, y, rlon, rlat
  REAL(r_sngl) :: wk(7)
  INTEGER      :: i, j, k, im, l, ll, ierr, iunit, iunit2, kmax
  INTEGER(8)   :: ncand, nkeep, nkeep_echo, ntrue_echo
  INTEGER(8), ALLOCATABLE :: nkeep_k(:)
  CHARACTER(5) :: cmem
  CHARACTER(10):: stdoutf = 'SYNO-00000'
  LOGICAL      :: ex
  REAL(r_size) :: rtimer00, rtimer

  CALL CPU_TIME(rtimer00)
  CALL initialize_mpi
  WRITE(stdoutf(6:10),'(I5.5)') myrank
  OPEN(6,FILE=stdoutf)

  INQUIRE(FILE='synth_obs.namelist',EXIST=ex)
  IF( ex )THEN
    OPEN(54,FILE='synth_obs.namelist')
    READ(54,NML=synth_obs,IOSTAT=ierr)
    CLOSE(54)
    IF( ierr > 0 )THEN
      WRITE(6,*)'ERROR: could not read &synth_obs in synth_obs.namelist'
      CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    ENDIF
  ENDIF
  WRITE(6,NML=synth_obs)

  !Same detection limit as the LETKF (set_letkf_obs, process_radar_obs).
  minrefdbz = dbz_min
  minref    = 10.0d0**( minrefdbz / 10.0d0 )
  use_pseudorh = .false.

  CALL set_common_wrf(TRIM(truth_file))
  WRITE(6,'(A,3I6)') 'Grid nlon, nlat, nlev: ', nlon, nlat, nlev

  ALLOCATE( v3d(nlon,nlat,nlev,nv3d), v2d(nlon,nlat,nv2d) )
  ALLOCATE( nact_loc(nlon,nlat,nlev) )
  nact_loc = 0
  elm  = REAL(id_reflectivity_obs,r_size)
  zero = 0.0d0

  !--------------------------------------------------------------------
  ! 1. Members with echo at every grid point
  !--------------------------------------------------------------------
  ll = CEILING( REAL(nens) / REAL(nprocs) )
  DO l = 1, ll
    im = myrank + 1 + (l-1)*nprocs
    IF( im > nens ) EXIT
    WRITE(cmem,'(I5.5)') im
    WRITE(6,'(A,I5,2A)') 'MYRANK ', myrank, ' reads ', TRIM(ens_prefix)//cmem
    CALL read_grd(TRIM(ens_prefix)//cmem,1,v3d,v2d)
    DO k = kmin, nlev-2
      DO j = 3, nlat-2
        DO i = 3, nlon-2
          typ  = radar_wavelength
          odat = 0.0d0
          CALL Trans_XtoY(elm,typ,zero,zero,odat,REAL(i,r_size),REAL(j,r_size),REAL(k,r_size), &
                          zero,zero,v3d,v2d,href)
          IF( href > minref ) nact_loc(i,j,k) = nact_loc(i,j,k) + 1
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  IF( myrank == 0 ) THEN
    ALLOCATE( nact(nlon,nlat,nlev) )
  ELSE
    ALLOCATE( nact(1,1,1) )
  ENDIF
  CALL MPI_REDUCE(nact_loc,nact,nlon*nlat*nlev,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  DEALLOCATE( nact_loc )

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(ENSEMBLE):', rtimer, rtimer-rtimer00
  rtimer00 = rtimer

  !--------------------------------------------------------------------
  ! 2-3. Truth, noise and output (rank 0)
  !--------------------------------------------------------------------
  IF( myrank == 0 ) THEN
    CALL read_grd(TRIM(truth_file),1,v3d,v2d)
    ALLOCATE( zt(nlon,nlat,nlev), ht(nlon,nlat,nlev), nkeep_k(nlev) )
    zt = v3d(:,:,:,iv3d_ph) / gg
    ht = dbz_min
    DO k = kmin, nlev-2
      DO j = 3, nlat-2
        DO i = 3, nlon-2
          typ  = radar_wavelength
          odat = 0.0d0
          CALL Trans_XtoY(elm,typ,zero,zero,odat,REAL(i,r_size),REAL(j,r_size),REAL(k,r_size), &
                          zero,zero,v3d,v2d,href)
          IF( href > minref ) ht(i,j,k) = 10.0d0*LOG10(href)
        ENDDO
      ENDDO
    ENDDO

    iunit  = 71
    iunit2 = 72
    OPEN(iunit, FILE=TRIM(obs_file),      FORM='unformatted',ACCESS='sequential',STATUS='replace')
    OPEN(iunit2,FILE=TRIM(truth_obs_file),FORM='unformatted',ACCESS='sequential',STATUS='replace')
    CALL ij_to_latlon(projection,REAL(nlon/2,r_size),REAL(nlat/2,r_size),rlat,rlon)
    WRITE(iunit)  REAL(rlon,r_sngl)
    WRITE(iunit)  REAL(rlat,r_sngl)
    WRITE(iunit)  0.0e0
    WRITE(iunit2) REAL(rlon,r_sngl)
    WRITE(iunit2) REAL(rlat,r_sngl)
    WRITE(iunit2) 0.0e0

    CALL init_gen_rand(seed)
    ncand = 0 ; nkeep = 0 ; nkeep_echo = 0 ; ntrue_echo = 0 ; nkeep_k = 0 ; kmax = 0
    DO k = kmin, nlev-2, skip_k
      DO j = 3, nlat-2, skip_ij
        DO i = 3, nlon-2, skip_ij
          IF( zt(i,j,k) > zmax ) CYCLE
          ncand = ncand + 1
          eps = gauss()                 !one draw per candidate, kept or not
          IF( ht(i,j,k) > dbz_min ) ntrue_echo = ntrue_echo + 1
          IF( nact(i,j,k) < min_active ) CYCLE
          y = MAX( ht(i,j,k) + obs_err*eps , dbz_min )
          CALL ij_to_latlon(projection,REAL(i,r_size),REAL(j,r_size),rlat,rlon)
          wk(1) = REAL(id_reflectivity_obs,r_sngl)
          wk(2) = REAL(rlon,r_sngl)
          wk(3) = REAL(rlat,r_sngl)
          wk(4) = REAL(zt(i,j,k),r_sngl)
          wk(5) = REAL(to_z(y),r_sngl)
          wk(6) = REAL(obs_err,r_sngl)
          wk(7) = REAL(radar_wavelength,r_sngl)
          WRITE(iunit) wk
          wk(5) = REAL(to_z(ht(i,j,k)),r_sngl)
          WRITE(iunit2) wk
          nkeep = nkeep + 1
          nkeep_k(k) = nkeep_k(k) + 1
          kmax = MAX(kmax,k)
          IF( ht(i,j,k) > dbz_min ) nkeep_echo = nkeep_echo + 1
        ENDDO
      ENDDO
    ENDDO
    CLOSE(iunit)
    CLOSE(iunit2)

    WRITE(6,'(A)')      '============================================='
    WRITE(6,'(A,I12)')  ' Candidate points            : ', ncand
    WRITE(6,'(A,I12)')  '   with echo in the truth    : ', ntrue_echo
    WRITE(6,'(A,I12)')  ' Written (n_active >= min)   : ', nkeep
    WRITE(6,'(A,I12)')  '   with echo in the truth    : ', nkeep_echo
    WRITE(6,'(A,I6)')   ' Highest level written       : ', kmax
    WRITE(6,'(A,F10.2)')' Approx. LETKF obshdxf [GB]  : ', REAL(nkeep,r_size)*REAL(nens,r_size)*8.0d0/1.0d9
    WRITE(6,'(A)')      ' Level, points written, mean z [m]'
    DO k = kmin, nlev-2
      IF( nkeep_k(k) > 0 ) WRITE(6,'(I6,I12,F10.1)') k, nkeep_k(k), SUM(zt(3:nlon-2,3:nlat-2,k))/REAL((nlon-4)*(nlat-4),r_size)
    ENDDO
    WRITE(6,'(A)')      '============================================='
  ENDIF

  CALL CPU_TIME(rtimer)
  WRITE(6,'(A,2F10.2)') '### TIMER(TRUTH+WRITE):', rtimer, rtimer-rtimer00
  WRITE(6,'(A)') ' SYNTH_OBS ends SUCCESFULLY !!!!'
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL finalize_mpi

CONTAINS

  !Box-Muller, one N(0,1) number per call, from the SFMT stream.
  REAL(r_size) FUNCTION gauss()
    REAL(r_dble) :: u1, u2
    u1 = genrand_res53()
    DO WHILE( u1 <= 0.0d0 )
      u1 = genrand_res53()
    ENDDO
    u2 = genrand_res53()
    gauss = SQRT( -2.0d0*LOG(u1) ) * COS( 2.0d0*pi*u2 )
  END FUNCTION gauss

  !dBZ -> Z [mm6 m-3]; values at or below the detection limit go to 0.
  REAL(r_size) FUNCTION to_z(dbz)
    REAL(r_size), INTENT(IN) :: dbz
    IF( dbz > dbz_min )THEN
      to_z = 10.0d0**( dbz / 10.0d0 )
    ELSE
      to_z = 0.0d0
    ENDIF
  END FUNCTION to_z

END PROGRAM make_synth_obs
