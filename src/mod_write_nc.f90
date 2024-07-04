MODULE mod_write_nc
  !!
  USE netcdf
  !!
  IMPLICIT none
  !!
  PUBLIC
  !!
  CHARACTER(len=200), PARAMETER :: &
       &    cv_lon = 'x',  &
       &    cv_lat = 'y',  &
       &    cv_t   = 'time'
  !
  !&    cv_lon = 'lon',  &
  !&    cv_lat = 'lat',  &
  !&    cv_t   = 'time'
  !!
  !!
  !!
  PUBLIC :: WRITE_NC
  !!
  !!
CONTAINS
  !!
  !!
  SUBROUTINE WRITE_NC(ni, nj, nt, vlon, vlat, vtime, xfld, cf_out, cv_out)
    !!
    INTEGER,                   INTENT(in) :: ni, nj, nt
    REAL, DIMENSION(ni),       INTENT(in) :: vlon
    REAL, DIMENSION(nj),       INTENT(in) :: vlat
    REAL, DIMENSION(nt),       INTENT(in) :: vtime
    REAL, DIMENSION(ni,nj,nt), INTENT(in) :: xfld
    CHARACTER(len=*),          INTENT(in) :: cf_out, cv_out
    !!
    !!
    !!
    INTEGER :: &
         &   ji, jj, jt, &
         &   ierr, id_out, id_lon, id_lat, id_t, id_x, id_y, id_time, id_field1
    !!
    !!
    !!
    !! Create file with ID id_out
    ierr = NF90_CREATE(cf_out, nf90_Write, ncid = id_out) ; call what_error(ierr)
    !! NF90_HDF5, nf90_noclobber, ....
    !!
    !! DIMENSIONS
    !! ~~~~~~~~~~
    !!
    !! Create longitude DIMENSION with ID id_x:
    ierr = NF90_DEF_DIM(id_out, trim(cv_lon), ni, id_x) ; call what_error(ierr)
    !!
    !! Create latitude DIMENSION with ID id_y:
    ierr = NF90_DEF_DIM(id_out, trim(cv_lat), nj, id_y) ; call what_error(ierr)
    !!
    !! Create time DIMENSION with ID id_t:
    ierr = NF90_DEF_DIM(id_out, trim(cv_t), nf90_unlimited, id_t) ; call what_error(ierr)
    !!
    !!
    !!
    !!
    !!
    !! VARIABLES
    !! ~~~~~~~~~~
    !!
    !! Create longitude VARIABLE with ID id_lon:
    ierr = NF90_DEF_VAR(id_out, trim(cv_lon), nf90_double, id_x, id_lon) ; call what_error(ierr)
    !!
    !! Create latitude VARIABLE with ID id_lat:
    ierr = NF90_DEF_VAR(id_out, trim(cv_lat), nf90_double, id_y, id_lat) ; call what_error(ierr)
    !!
    !! Create time VARIABLE with ID id_time:
    ierr = NF90_DEF_VAR(id_out, trim(cv_t),   nf90_double, id_t, id_time) ; call what_error(ierr)
    !!
    !!
    !! Create OUTPUT FIELD VARIABLE (3D -> ni,nj,nt):
    ierr = NF90_DEF_VAR(id_out, trim(cv_out), nf90_double, (/id_x, id_y, id_t /), id_field1)
    call what_error(ierr)
    !!
    !!
    !! ATRIBUTES, for each variable and global file
    !!
    !! => NF90_PUT_ATT(...)
    !! .......................
    !! .......................
    !! .......................
    !!
    !!
    !! DEFINITION IS OVER
    !! ~~~~~~~~~~~~~~~~~~
    ierr = NF90_ENDDEF(id_out) ; call what_error(ierr)
    !!
    !!
    !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !!
    !!
    !!
    !! TIME TO FILL THE VARIABLES
    !! ~~~~~~~~~~~~~~~~~~~~~~~~~~
    !!
    !! Writing longitude:
    ierr = NF90_PUT_VAR(id_out, id_lon, VLON) ; call what_error(ierr)
    !!
    !! Writing latitude:
    ierr = NF90_PUT_VAR(id_out, id_lat, VLAT) ; call what_error(ierr)
    !!
    !! Writing latitude:
    ierr = NF90_PUT_VAR(id_out, id_time, vtime) ; call what_error(ierr)
    !!
    !!
    DO jt = 1, nt
       !!
       !ierr = NF90_PUT_VAR(id_out, id_time, vtime(jt), start = jt, count = 1 )
       !call what_error(ierr)
       !!
       ierr = NF90_PUT_VAR(id_out, id_field1, xfld(:,:,jt), start = (/ 1, 1, jt /), count = (/ ni, nj, 1 /) )
       call what_error(ierr)
       !!
    END DO
    !!
    !!
    ierr = NF90_CLOSE(id_out) ; call what_error(ierr)
    !!
    !!
  END SUBROUTINE WRITE_NC
  !!
  !!
  !!
  !!
  !!
  !!
  SUBROUTINE what_error(ierror)
    !!
    IMPLICIT none
    !!
    INTEGER, INTENT(in) :: ierror
    !!
    IF ( ierror /= 0 ) THEN
       PRINT *, 'There was an error, # = ', ierror
       STOP
    END IF
    !!
  END SUBROUTINE what_error
  !!
  !!
  !!
END MODULE mod_write_nc
