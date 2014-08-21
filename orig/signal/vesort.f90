PROGRAM vesort

! Sort vector read from some file with many columns
! Then compute median and quartiles using sorted vector

REAL, DIMENSION(1000000) :: v
REAL :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20
REAL :: enp,vmedian,vquart1,vquart2
INTEGER :: i,ico,np,nenp,menp

CHARACTER (LEN=132) :: cline
CHARACTER (LEN=50) :: ifile

CALL GETARG(1,cline)
READ(cline,*) ifile
CALL GETARG(2,cline)
READ(cline,*) ico

OPEN(10,FILE=ifile)
OPEN(11,FILE='vesort.dat')

np = 0
DO i = 1, 1000000
  READ(10,*,END=10) v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20
  IF (ico == 1) v(i) = v1
  IF (ico == 2) v(i) = v2
  IF (ico == 3) v(i) = v3
  IF (ico == 4) v(i) = v4
  IF (ico == 5) v(i) = v5
  IF (ico == 6) v(i) = v6
  IF (ico == 7) v(i) = v7
  IF (ico == 8) v(i) = v8
  IF (ico == 9) v(i) = v9
  IF (ico == 10) v(i) = v10
  IF (ico == 11) v(i) = v11
  IF (ico == 12) v(i) = v12
  IF (ico == 13) v(i) = v13
  IF (ico == 14) v(i) = v14
  IF (ico == 15) v(i) = v15
  IF (ico == 16) v(i) = v16
  IF (ico == 17) v(i) = v17
!  WRITE(*,*) i,v(i)
  np = np + 1
END DO
10 CONTINUE
CLOSE(10)

WRITE(*,*) 'Data points:', np

vmin = 1000000.0
vmax = -1000000.0
imin = 0
imax = 0
DO i = 1, np
  IF (v(i) .LT. vmin) THEN
    vmin = v(i)
    imin = i
  END IF
  IF (v(i) .GT. vmax) THEN
    vmax = v(i)
    imax = i
  END IF
END DO
WRITE(*,*) 'Vector extrema:', imin,vmin,imax,vmax

CALL tisort(v,np)

DO i = 1, np
  WRITE(11,'(I9,1X,F9.3)') i, v(i)
END DO
CLOSE(11)

! Median and quartiles
enp = REAL(np)/2.0
nenp = ANINT(enp+0.1)
menp = ANINT(enp-0.1)
!WRITE(*,*) enp,nenp,menp
IF (nenp == menp) vmedian = 0.5*(v(nenp)+v(nenp+1))
IF (nenp > menp) vmedian = v(nenp)
enp = REAL(np)/4.0
nenp = ANINT(enp+0.1)
menp = ANINT(enp-0.1)
!WRITE(*,*) enp,nenp,menp
IF (nenp == menp) vquart1 = 0.5*(v(nenp)+v(nenp+1))
IF (nenp > menp) vquart1 = v(nenp)
enp = 0.75*REAL(np)
nenp = ANINT(enp+0.1)
menp = ANINT(enp-0.1)
!WRITE(*,*) enp,nenp,menp
IF (nenp == menp) vquart2 = 0.5*(v(nenp)+v(nenp+1))
IF (nenp > menp) vquart2 = v(nenp)
WRITE(*,*) 'Median', vmedian,vquart1,vquart2
  

CONTAINS

! Sort array in increasing order
SUBROUTINE tisort(art,na)

REAL, DIMENSION(100000), INTENT(INOUT) :: art
INTEGER, INTENT(IN) :: na
REAL :: tite

INTEGER :: i, imi(1), ime

DO i = 1, na-1
  imi = MINLOC(art(i:na))
  ime = imi(1)+i-1
  IF (ime .NE. i) THEN
    tite = art(i)
    art(i) = art(ime)
    art(ime) = tite
  END IF
END DO

RETURN

END SUBROUTINE tisort

END PROGRAM vesort
