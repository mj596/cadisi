 PROGRAM flasim

! Version of 25 September 2012
! Simulates/computes FlashCam signal reconstruction
! Initial version with photomultiplier as a photon detector

IMPLICIT NONE

REAL, DIMENSION(100000) :: t1,t2,tt,te,pude,hit,pufi,ta,fadi,ti,fadan,fadepi
REAL, DIMENSION(100000) :: fada,tu,fadu,fade,fades,fadp,fadep,fadesa,fadar,fadeco
REAL, DIMENSION(100000) :: fady,fasa,fast,fadc,fadin,wnoise,fadn,fasin,fadeni,fadb
REAL, DIMENSION(100000) :: nsbart,cerart,atime,ptime,spex,sper,am,fader,fadap,fadapi
REAL :: atimemin,atimemax,sumati,artime,detime,sitime,time0,uh,th,sh,ah,ytest
REAL :: desum,amsum,ares,deta,tsam,tstart,tadc,dt,onepe,tad,fasum,ai,tic,deti
REAL :: awe,starti,sisum,det1,det2,bi,detimin,tam1,tam2,dti,swe,resum,sampli
REAL :: recoti,sampsum,repe,refepe,pudemax,pufisum,pufimax,sampmax,isamsum,fisum
REAL :: isammax,misasum,misamax,tave,dth,af1,af2,af3,wnom,wnos,wnomean,wnoisa
REAL :: tno,nsbrate,tstop,tra,tde,fwhm1,fwhm2,fwhmsum,fwhmtime,refeped,refepes
REAL :: refepedsum,repede,bgsumde,bgsumsh,spet,samplin,recotin,recota,recotu
REAL :: decotin,samplide,posisum,nesisum,pomean,nemean,anesisum,anemean
REAL :: dif1,dif2,dif3,dif4,fasumi,defa,signorm

INTEGER, DIMENSION(100000) :: nf,nno
INTEGER, DIMENSION(1000) :: me
INTEGER :: i,k,j,nc,istart,istop,ne,nseed,j0,j1,j2,j3,jj1,jj2,jj3,na,ii,ni,nt
INTEGER :: nsbnum,np,ihm1,ihm2,nbgde,nbgsh,nspe,npos,nneg,ncer
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

REAL :: cerate,nsrate,ton,toff,events,photons,pi,sigma1p,gno,gen,gex,celec1p
REAL :: sigma0,celec0,gex0,gno0,gen0,darate,sigma1,celec1,gnop,gexp,genp,sigmap
REAL :: alpha,we,ala,beta,era,alpa,celecn,sigman,ar1,ar2,sia,tevents,wep
REAL :: sump,sumpe,mumax,pmax,mumean,cemu,sumap,mumedian,sigmamu,clower,cupper
REAL :: sigmal,sigmau,sumpi,sumpa,mimean,sigmami,mimedian,nsmean,nssigma
REAL :: cnsb,wnsb,npnsb,sumns,sumfadis,ranu,rana,maxfadis,den0,dem0
REAL :: sigma02,sigma12,telaps1,telaps2
REAL, DIMENSION(1000) :: nfoto,pnsb
REAL, DIMENSION(20001) :: p,px,pn,pxn,wo,rolo,bemax,rpc,rpb,fadis,plo,plox
REAL, DIMENSION(20) :: runtim1,runtim2
INTEGER :: non,noff,kk,mu,kkmax,kb,pshift,nmu,radc,clock
INTEGER :: fadc1,fadc2,ietime
CHARACTER (LEN=132) :: cline
CHARACTER (LEN=10) :: disfile,resfile
CHARACTER (LEN=7) :: sodir
CHARACTER (LEN=1) :: cleft,cun

call ETIME(runtim1,telaps1)

disfile = 'spesim.dat'
resfile = 'spesim.res'
cleft = '/'
cun = '_'

pi = 4.0*ATAN(1.0)
clower = (1.0-0.6827)/2.0
cupper = 1.0 - clower

! Parameters
! Cherenkov photons
ncer = 10
! Mean arrival time of Cherenkov photons
time0 = 100.0
! Sampling period
tsam = 4.0
!tsam = 1.0
nseed = 17
! PMT amplitude corresponding to 1 photoelectron
! 250 MHz
onepe = 28.418
! 1 GHz
!onepe = 14.312
! White noise parameters
wnom = 0.0
!wnos = 4.7
!wnos = 20.0
!wnos = 0.001
! Cracow: 0.15 mV noise for 1 mV of one p.e. signal
wnos = 4.263
! NSB rate MHz
!nsbrate = 250.0
!nsbrate = 125.0
! Thermal noise of GPAD, 6 MHz
nsbrate = 6.0
! Derivative added to the averaged signal with factor defa
!defa = 0.333333
! PMT 250 Mhz
defa = 0.9
! PMT 1000 MHz
!defa = 0.63
! Single photon signal normalization
! PMT 250 MHz
!signorm = 6.73
! PMT 1000 MHz
!signorm = 2.70
! GAPD
signorm = 19.7789

CALL RANDOM_SEED(size = nseed)
ALLOCATE(seed(nseed))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i, i = 1, nseed) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)

! Read photon bunch arrival times
! Set reference time (100 ns) to the weighted mean time of the bunch
! Reference time: photons, not ADC
tstart = 0.0
tstop = 1000.0
dt = 0.1
dti = 0.01
nt = NINT((tstop-tstart)/dt)
!ni = NINT((tstop-tstart)/dti)

!! WRITE(*,*) 'Coarse and fine bins:', nt, ni

OPEN(10,FILE='flasim.test')
OPEN(11,FILE='white.noise')
! Single photoelectron spectrum with cross-talk and afterpulses
! fitted to the mppc_uknown.dat data
OPEN(12,FILE='mppcunk.dat')
!OPEN(13,FILE='speresp.test')
OPEN(14,FILE='photon.bunches')
OPEN(19,FILE='flasim.results')
OPEN(20,FILE='flasim.dist')
!OPEN(11,FILE='photons10.dat')
OPEN(21,FILE='cherenkov.photons')
OPEN(22,FILE='nsb.photons')
OPEN(23,FILE='detector.signals')
OPEN(24,FILE='detector.pulse')
OPEN(25,FILE='amplifier.pulse')
OPEN(26,FILE='fadc.samples')
OPEN(27,FILE='fadc.inter')
OPEN(28,FILE='fadc.aver')
OPEN(29,FILE='fadc.averin')
OPEN(30,FILE='fadc.dave')
OPEN(31,FILE='fadc.daves')
OPEN(32,FILE='fadc.daver')
OPEN(33,FILE='fadc.davesa')
OPEN(34,FILE='fadc.deco')
OPEN(35,FILE='fadc.shaped')
OPEN(36,FILE='fadc.shapin')

nspe = 0
DO i = 1, 10000
  READ(12,*,END=3) spex(i),sper(i)
  nspe = nspe + 1
END DO
3 CONTINUE
CLOSE(12)

!! WRITE(*,*) 'SPE file lines:', nspe

! Test if speresp function works properly
!DO i = 1, 50000
!  spet = speresp(1.0,spex,sper,nspe)
!  WRITE(13,*) spet
!END DO
!CLOSE(13)

! Cherenkov photons arriveal time read from the sim_telarray log
! or selected according to a box distribution (uniform with gaussian tails)

!DO i = 1, 1000
! photons10.dat
!  READ(11,*,END=10) atime
!END DO
!10 CONTINUE
!CLOSE(11)

!DO i = 1, 1000
!  hit(i) = 0.0
!END DO

!DO i = 1, 1000000
!  artime = rano(100.0,0.5)
!  DO j = 1, 1000
!    t1(j) = 95.0 + REAL(j)*0.01
!    IF (artime > t1(j) .AND. artime < t1(j)+0.01) THEN
!      hit(j) = hit(j) + 1.0
!    END IF
!  END DO
!END DO

!DO i = 1, 1000
!  hit(i) = hit(i)/30.0
!  WRITE(10,*) t1(i),hit(i)
!END DO
!CLOSE(10)

! Generate random arrival times of photons
! Add jitter time
DO i = 1, ncer
!  atime(i) = bodi(100.0,1.4)
  atime(i) = bodi(100.0,10.0)
!  WRITE(*,*) atime(i)
  atime(i) = atime(i) + rano(0.0,0.5)
!!  WRITE(*,*) atime(i)
END DO 

! Instrumental resolution tests, all photons arrive at t = 100 ns
!DO i = 1, ncer
!  atime(i) = 100.0
!END DO

atimemin = 1000.0
atimemax = 0.0
DO i = 1, ncer
  IF (atime(i) < atimemin) atimemin = atime(i) 
  IF (atime(i) > atimemax) atimemax = atime(i)
END DO
!WRITE(*,*) 'Time range:', atimemin,atimemax

nc = 0
sumati = 0.0
DO i = 1, 1000
  IF (atime(i) > 0.0) THEN
    atime(i) = atime(i) - atimemin
    sumati = sumati + atime(i)
!    WRITE(*,*) i, atime(i)
    nc = nc + 1
  END IF
END DO

artime = sumati/REAL(nc)
detime = atimemax - atimemin

sumati = 0.0
DO i = 1, nc
  sumati = sumati + (atime(i)-artime)**2.0
END DO

sitime = SQRT(sumati/REAL(nc))

WRITE(19,*) 'Photons:', nc, artime, sitime, detime

DO i = 1, nc
  cerart(i) = time0 + atime(i) - artime
END DO

DO i = 1, nc
! cherenkov.photons
  WRITE(21,*) cerart(i)
END DO
CLOSE(21)

! Simulate NSB photons arrival times

CALL nsbat(nsbrate,tstart,tstop,nsbnum,nsbart)

CALL tisort(nsbart,nsbnum)

!! WRITE(*,*) 'NSB photons:', nsbnum

!nsbnum = 0
!DO i = 1, 10000
!  nsbnum = nsbnum + 1
!  READ(22,*,END=10) nsbart(i)
!END DO
!10 CONTINUE

DO i = 1, nsbnum
! nsb.photons
  WRITE(22,*) nsbart(i)
!!  WRITE(*,*) i, nsbart(i)
END DO
CLOSE(22)

! Merge Cherenkov and NSB photons into photon list
! DIFFERENT QUANTUM EFFICIENCY (DUE TO DIFFERENT SPECTRA) CAN BE TAKEN INTO ACCOUNT
np = 1
DO i = 1, nsbnum
  ptime(np) = nsbart(i)
  np = np + 1
END DO

DO i = 1, nc
  ptime(np) = cerart(i)
  np = np + 1
END DO
np = np - 1

CALL tisort(ptime,np)

!! WRITE(*,*) 'Total photons:', np
!! DO i = 1, np
!!   WRITE(*,*) i, ptime(i)
!! END DO

! Investigated time period: 70--170 ns

!nt = 2000
DO i = 1, nt
  t1(i) = tstart + REAL(i-1)*dt
  t2(i) = t1(i) + dt
  tt(i) = t1(i) + 0.5*dt
  nf(i) = 0
  pude(i) = 0.0
END DO

! Sum photons over time bins

DO i = 1, nt
  DO j = 1, np
    IF (ptime(j) >= t1(i) .AND. ptime(j) < t2(i)) nf(i) = nf(i) + 1
  END DO
END DO

ne = 0
!istart = 1
!istop = 1
DO i = 1, nt
!  IF (ne < 1) istart = i
  IF (nf(i) > 0) THEN
!    WRITE(*,*) i, tt(i), nf(i)
! detector.signals
    WRITE(23,*) i, tt(i), nf(i)
    ne = ne + 1
    te(ne) = tt(i)
    me(ne) = nf(i)
!    istop = i
  END IF
END DO
CLOSE(23)
!istart = istart+1
! Add 20 ns for PMT signal
!istop = istop + 200
! Add 150 ns for GAPD signal
!istop = istop + 1500
!IF (istop > nt) istop = nt
!WRITE(*,*) 'istart,istop', istart,istop

! Signal amplitude blurred according to single photoelectron spectrum
! Each photon amplitude blurred separately and then summed up
DO j = 1, ne
  am(j) = 0.0
  DO k = 1, me(j)
    am(j) = am(j) + speresp(onepe,spex,sper,nspe)
!    am(j) = REAL(me(j))*onepe
  END DO
! photon.bunches
  WRITE(14,*) te(j),me(j)
END DO
CLOSE(14)

! Sum elementary pulses into detector pulse ('signal' in figure on page 37)

! 250 MHz
sh = 1.4154
uh = 0.83725
dth = 2.95
! Silicon detectors pulse shape
tra = 1.12
tde = 25.33
dth = 0.0
! 1 GHz
!sh = 1.425
!uh = 0.90954
!dth = 3.01
!DO i = istart, istop
DO i = 1, nt
  DO j = 1, ne
    IF (te(j) <= tt(i)) THEN
      th = te(j) + dth
      ah = REAL(me(j))*onepe
! PMT signal
!      pude(i) = pude(i) + hyga(tt(i),ah,th,sh,uh)
! GAPD signal
      pude(i) = pude(i) + dexe(tt(i),am(j),th,tra,tde)
    END IF
  END DO
END DO

desum = 0.0
DO i = 1, nt
  desum = desum + pude(i)
!  IF (tt(i) < 100.0) ytest = 0.0
!  IF (tt(i) >= 100.0) ytest = hyga(tt(i),284.18,102.95,sh,uh)
! detector.pulse
!  WRITE(24,*) tt(i), pude(i)
!  WRITE(10,*) tt(i), ytest
END DO
!CLOSE(24)
!CLOSE(10)

! Reference pulse integral
resum = desum*dt
pudemax = MAXVAL(pude)

WRITE(19,*) 'Detector pulse integral and maximum:', resum, pudemax
!! WRITE(*,*) 'Detector pulse integral and maximum:', resum, pudemax

! Photons arrive to ADC at random ADC phase
! Sampling starts at tstart + 0.75*tsam (see sampling algorithm)
! Therefore ADC start time is in the range (0.75*tsam,1.75*tsam)
tadc = bodi(tstart+1.25*tsam,tsam)
!tadc = tstart
!tadc = tstart+0.5*tsam
!tadc = tstart+tsam
!tadc = 71.200005
!tadc = 71.3
!tadc = tstart + 0.9*tsam
! Basic phase (for signal model)
!tadc = tstart + 0.82*tsam
na = 1
tad = tadc
DO WHILE (tad < tt(nt)-0.5*tsam)
  ta(na) = tadc + tsam*(REAL(na)-0.5)
  tad = ta(na)
  na = na + 1
END DO
na = na - 1
!! WRITE(*,*) 'ADC start time, na:', tadc, na, tad

! Fine interpolation bins
ni = 1
tic = tadc
DO WHILE (tic < tad)
  ti(ni) = tadc + REAL(ni)*dti
  tic = ti(ni)
  ni = ni + 1
END DO
ni = ni-1
!! WRITE(*,*) 'Interpolation bins ni:', ni,ti(1),tic

! White noise, normal distribution of amplitudes
amsum = 0.0
fasum = 0.0
DO i = 1, 1000000
  tno = tstart + REAL(i)*0.001
  wnoisa = rano(wnom,wnos)
!  amsum = amsum + 0.0001*ABS(wnoisa)
!  fasum = fasum + 0.0001
  amsum = amsum + ABS(wnoisa)*ABS(wnoisa)
  DO j = 1, nt
    IF (tno >= tt(j)-0.5*dt .AND. tno < tt(j)+0.5*dt) THEN
! Noise normalization, fine bin 0.0001 ns divided by dt = 0.1 ns
      wnoise(j) = wnoise(j) + 0.01*wnoisa
!      nno(j) = nno(j) + 1
    END IF
  END DO
END DO
! White noise maplitude per ns
!wnomean = amsum/fasum
!wnomean = amsum/(1000000*0.0001)
wnomean = SQRT(amsum/1000000)
WRITE(*,*) 'White noise integral and amplitude (0.0001 ns):', amsum, wnomean

! white.noise added to detector signal
!amsum = 0.0
!fasum = 0.0
DO i = 1, nt
!  wnoise(i) = wnoise(i)/REAL(nno(i))
!  WRITE(11,*) tt(i),wnoise(i)
! No noise
!  pude(i) = pude(i) + onepe*wnoise(i)
!      IF (i == 310) WRITE(*,*) 'pude', i,tt(i),pude(i)
!  amsum = amsum + dt*ABS(wnoise(i))
!  fasum = fasum + dt
! detector.pulse
  WRITE(24,*) tt(i), pude(i)
!  WRITE(10,*) tt(i), ytest
END DO
CLOSE(24)
!CLOSE(11)
!wnomean = onepe*amsum/fasum
!WRITE(*,*) 'White noise integral and amplitude (0.1 ns):', amsum, wnomean

! Convolve the detector pulse with the amplifier/filters response

! 250 MHz
af1 = 6.5
af2 = 1.1
af3 = 4.5
! 1 GHz
!af1 = 3.0
!af2 = 10.0
!af3 = 1.15
fisum = 0.0
DO i = 1, nt
  IF (i == 1) pufi(i) = 0.0
  IF (i > 1) pufi(i) = amfire(i,pude,tt,af1,af2,af3)
  fisum = fisum + pufi(i)
  deta = tt(i)-tt(1)
  ares = 0.0
  ares = 1*(af1+af2*deta)*exp(-deta/af3)
  WRITE(10,*) deta, ares
END DO

! white.noise added to amplified and filtered signal
amsum = 0.0
fasum = 0.0
DO i = 1, nt
!  wnoise(i) = wnoise(i)/REAL(nno(i))
  WRITE(11,*) tt(i),wnoise(i)
! No noise
  pufi(i) = pufi(i) + onepe*wnoise(i)
!      IF (i == 310) WRITE(*,*) 'pude', i,tt(i),pude(i)
  amsum = amsum + dt*ABS(wnoise(i))
  fasum = fasum + dt
END DO
CLOSE(11)
wnomean = onepe*amsum/fasum
WRITE(*,*) 'White noise integral and amplitude (0.1 ns):', amsum, wnomean

DO i = 1, nt
  pufi(i) = desum*pufi(i)/fisum
! amplifier.pulse
  WRITE(25,*) tt(i), pufi(i)
END DO
CLOSE(25)
CLOSE(10)
pufisum = SUM(pufi)*dt
pufimax = MAXVAL(pufi)

WRITE(19,*) 'Amplifier pulse integral and maximum:', pufisum, pufimax

! Amplitude and time at signal maximum
sampli = -100000.0
DO i = 1, nt
  IF (pufi(i) > sampli) THEN
    sampli = pufi(i)
    recota = tt(i)
  END IF
END DO

WRITE(19,*) 'Amplifier signal, arrival time:', recota

!! WRITE(*,*) 'Amplifier pulse integral and maximum:', pufisum, pufimax

! Reference signal amplitude for tadc = tstart + 0.5*tsam (maximal signal)
! Should be the mean of distibution of results
! Deconvolved signal 
!refeped = 12.353
refeped = 43.43242
refepedsum = 86.972
! Shaped signal
!refepes = 7.429
refepes = 30.3221

! FADC sampling (SAMPLING, NOT REBINNING, SOME INFORMATION LOST)
! Average over a half clock cycle, starting 0.75 cycle before ta(i)
! because average is done over past pulse (-0.25,+0.25)
! at the ADC time bin start
! Rounding to the nearest integer
DO i = 1, na
  j1 = ANINT((ta(i)-0.75*tsam-tstart)/dt)+1
!  j1 = ii - ANINT(0.25*tsam/dt)
  IF (j1 < 1) j1 = 1
  IF (j1 > nt) j1 = nt
  j2 = ANINT((ta(i)-0.25*tsam-tstart)/dt)
  IF (j2 > nt) j2 = nt
  tam1 = ta(i)-0.75*tsam
  tam2 = ta(i)-0.25*tsam
!!  IF (i == 1)  WRITE(*,*) 'Sampling:',i,j1,j2,ta(i),tam1,tam2,tt(j1),tt(j2)
  amsum = 0.0
  DO j=j1,j2
    amsum = amsum + pufi(j)
!    IF (i > 8 .AND. i < 11) WRITE(*,*) j,tt(j),pufi(j)
  END DO
  fasum = amsum/REAL(j2-j1+1)
  fadc(i) = REAL(ANINT(fasum))
!  fadc(i) = fasum
! fadc.samples
  WRITE(26,*) ta(i),fadc(i)  
END DO
CLOSE(26)

! FADC samples time is SHIFTED AHEAD by 0.5 tsam
! Transform to the real time for further signal analysis
sampsum = 0.0
DO i = 1, na
  ta(i) = ta(i) - 0.5*tsam
  sampsum = sampsum + fadc(i)*tsam
END DO
sampmax = MAXVAL(fadc)

WRITE(19,*) 'FADC samples integral and maximum:', sampsum, sampmax
!! WRITE(*,*) 'FADC samples integral and maximum:', sampsum, sampmax

! Smoothing the ADC counts with a moving average of 1 clock cycle

! First linear interpolation is done to obtain ADC amplitude for 1 ns bins
! Linear and parabolic interpolation of the sampled signal

!CALL linin(na,ni,ta,ti,fadc,fadi)
!CALL parin(na,ni,ta,ti,fadc,fadp)
CALL herin(na,ni,ta,ti,fadc,fadi)
!CALL curin(na,ni,ta,ti,fadc,fadp)
!CALL porin(na,ni,ta,ti,fadc,fadp)

DO i = 1, ni
!  fadapi(i) = fadp(i)
! fadc.inter
  WRITE(27,*) ti(i),fadi(i),fadi(i)
END DO
CLOSE(27)  

isamsum = SUM(fadi)*dti
isammax = MAXVAL(fadi)

WRITE(19,*) 'Linear interpolation of FADC samples integral and maximum:', isamsum, isammax

! Moving average of interpolated ADC signal

!CALL movave(ni,nt,ti,tt,fadi,fada,1.0*tsam)
CALL movave(ni,nt,ti,tt,fadi,fada,1.0*tsam)
!CALL linin(nt,ni,tt,ti,fada,fadi)
!CALL movave(ni,nt,ti,tt,fadp,fadap,tsam)
!CALL linin(nt,ni,tt,ti,fada,fadi)
!CALL movave(ni,nt,ti,tt,fadi,fada,tsam)

fasum = 0.0
amsum = 0.0
DO i = 1, nt
! fadc.aver
  WRITE(28,*) tt(i),fada(i),fada(i)
!  IF (tt(i) > 200.0 .AND. tt(i) < 230.0) THEN
!    fasum = fasum + fada(i-1)-fada(i)
!    amsum = amsum + 1.0
!  END IF
END DO
CLOSE(28)
sampli = fasum/amsum
!WRITE(*,*) 'Mean fada decrease:', sampli,amsum

misasum = SUM(fada)*dt
misamax = MAXVAL(fada)

WRITE(19,*) 'Moving average of interpolated FADC samples integral and maximum:', misasum, misamax
 
! Smoothing the averaged FADC signal with another interpolation

DO i = 1, nt
  tu(i) = tt(i) - 0.5*dt
END DO

CALL herin(nt,ni,tt,ti,fada,fadi)
CALL herin(nt,nt,tt,tu,fada,fadu)
!CALL parin(nt,nt,tt,tu,fada,fadp)

DO i = 1, nt
! fadc.averin
  WRITE(29,*) tu(i),fadu(i),fadu(i)
END DO
CLOSE(29)      

! Amplitude and time at signal maximum
sampli = -100000.0
DO i = 1, ni
  IF (fadi(i) > sampli) THEN
    sampli = fadi(i)
    recotu = ti(i)
  END IF
END DO

WRITE(19,*) 'Averin signal, arrival time:', recotu

!WRITE(*,*) 'Signal integral:', sisum

! Simple derivative of the smoothed averaged FADC signal
amsum = 0.0
fasum = 0.0
DO i = 1, nt
!DO i = 1, ni
  fade(i) = 0.0
!  fadep(i) = 0.0
!  IF (i < nt) fade(i) = (fada(i+1)-fada(i-1))/dt
!  IF (i < nt) fadep(i) = (fady(i+1)-fady(i))/dt
!  IF (i < 1000) fade(i) = (fadp(i+1)-fadp(i))/dt
   IF (i > 2 .AND. i < nt-1) THEN
     fade(i) = (8.0*fada(i+1)+fada(i-2)-8.0*fada(i-1)-fada(i+2))/(12.0*dt)
   END IF
!   IF (i > 2 .AND. i < ni-1) THEN
!     fade(i) = (8.0*fadi(i+1)+fadi(i-2)-8.0*fadi(i-1)-fadi(i+2))/(12.0*dti)
!   END IF
!   IF (i > 2 .AND. i < nt-1) THEN
!     fadep(i) = (8.0*fadap(i+1)+fadap(i-2)-8.0*fadap(i-1)-fadap(i+2))/(12.0*dt)
!     fadep(i) = (fadap(i+1)-fadap(i))/dt
!     IF (tt(i) > 102.5 .AND. tt(i) < 103.5) THEN
!       dif1 = fadap(i-1)-fadap(i-2)
!       dif2 = fadap(i)-fadap(i-1)
!       dif3 = fadap(i+1)-fadap(i)
!       dif4 = fadap(i+2)-fadap(i+1)
!       WRITE(*,*) 'dave:',dif1,dif2,dif3,dif4,fadep(i)
!     END IF
!   END IF
!  WRITE(*,*) i, tu(i), tt(i), fadu(i), fade(i)
  amsum = amsum + dt*ABS(fade(i))
!  amsum = amsum + dti*ABS(fade(i))
!  fasum = fasum + dt*ABS(fadep(i))
END DO

WRITE(*,*) 'Derivative integral:', amsum,fasum
     
DO i = 1, nt
!DO i = 1, ni
!  fade(i) = resum*fade(i)/amsum
! Normalization from a single photon signal
fade(i) = signorm*fade(i)
!  fadep(i) = resum*fadep(i)/fasum
! fadc.dave
  WRITE(30,*) tt(i),fade(i),fade(i)
!  WRITE(30,*) ti(i),fade(i),fade(i)
END DO
CLOSE(30)      

! Fine interpolation of the averaged signal 

CALL herin(nt,ni,tt,ti,fade,fadi)
!CALL herin(ni,ni,ti,ti,fade,fadi)
!CALL linin(nt,ni,tu,ti,fadep,fadp)

!DO i = 1, ni
!  fadepi(i) = resum*fadepi(i)/fasumi
! flasim.dist
!  WRITE(20,*) ti(i),fadepi(i)
!END DO
!CLOSE(20)

! Moving average of interpolated ADC signal

! PMT signal, 250 MHz
!CALL movave(ni,nt,ti,tt,fadi,fadar,tsam)
CALL movave(ni,nt,ti,tt,fadi,fadar,1.0*tsam)
! GAPD signal
!CALL movave(ni,nt,ti,tt,fadi,fadar,1.0*tsam)
!CALL movave(ni,nt,ti,tt,fadp,fadan,2.0*tsam)

DO i = 1, nt
! fadc.daver
  WRITE(32,*) tt(i),fadar(i),fadar(i)
END DO
CLOSE(32)

! Interpolation of the averaged signal 

CALL herin(nt,ni,tu,ti,fadi,fades)

DO i = 1, nt
! fadc.daves
  WRITE(31,*) tu(i),fades(i)
END DO
CLOSE(31)      

! Moving average of interpolated ADC signal

CALL movave(nt,nt,tt,tu,fades,fadesa,0.5)

Do i = 1, nt
! fadc.davesa
  WRITE(33,*) tt(i),fadesa(i),j2-j1
END DO
CLOSE(33)

! Deconvolved signal: sum of averaged signal and averaged derivative
DO i = 1, nt
  fadeco(i) = 0.0
!  IF (i > 1) fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + fadesa(i)
! PMT signal
!  IF (i > 1) fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + fadar(i)
!  IF (i > 1) fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + fade(i)
! GAPD signal
  IF (i > 1) fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + defa*fadar(i)
!  IF (i > 1) fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + 0.5*fade(i)
! fadc.deco
  WRITE(34,*) tt(i),fadeco(i)
END DO
CLOSE(34)

! Positive signal level
posisum = 0.0
npos = 0
DO i = 1, nt
!  IF (tt(i) > 10.0 .AND. tt (i) < 80.0) THEN
!    posisum = posisum + fadu(i)
!    npos = npos + 1
!  END IF
  IF (tt(i) > 0.0 .AND. tt (i) < 900.0) THEN
    posisum = posisum + fadu(i)
    npos = npos + 1
  END IF
END DO

! Negative signal level
nesisum = 0.0
anesisum = 0.0
nneg = 0
DO i = 1, nt
!  IF (tt(i) > 10.0 .AND. tt (i) < 80.0) THEN
!    IF (fadar(i) <0.0) nesisum = nesisum + fadar(i)
!    IF (fadar(i) >= 0.0) anesisum = anesisum + ABS(fadar(i))
!    nneg = nneg + 1
!  END IF
  IF (tt(i) > 0.0 .AND. tt (i) < 900.0) THEN
    IF (fadar(i) < 0.0) nesisum = nesisum + fadar(i)
    IF (fadar(i) >= 0.0) anesisum = anesisum + ABS(fadar(i))
    nneg = nneg + 1
  END IF
END DO

pomean = posisum/REAL(npos)
nemean = nesisum/REAL(nneg)
anemean = anesisum/REAL(nneg)

WRITE(*,*) 'Positive and negative signal level:', pomean,nemean,anemean

! Amplitude and time at signal maximum
sampli = -100000.0
DO i = 1, nt
  IF (tt(i) > 420.0) EXIT
  IF (fadeco(i) > sampli) THEN
    sampli = fadeco(i)
    recoti = tt(i)
  END IF
END DO

! Fine interpolation for a better time resolution
CALL herin(nt,ni,tt,ti,fadeco,fadeni)

samplide = -100000.0
DO i = 1, ni
  IF (ti(i) > 420.0) EXIT
  IF (fadeni(i) > samplide) THEN
    samplide = fadeni(i)
    decotin = ti(i)
  END IF
END DO

! Amplitude and time from the (-0.5FWHM,+0.5FWHM) range
fwhm1 = 0.0
fwhm2 = 0.0
DO i = 1, nt
  IF (fadeco(i) >= 0.5*sampli .AND. fadeco(i-1) < 0.5*sampli) THEN
    fwhm1 = tt(i)
    ihm1 = (i)
  END IF
  IF (fadeco(i) < 0.5*sampli .AND. fadeco(i-1) >= 0.5*sampli) THEN
    fwhm2 = tt(i)
    ihm2 = i
  END IF
END DO

fwhmsum = 0.0
fwhmtime = 0.0
DO i = ihm1, ihm2
  fwhmsum = fwhmsum + fadeco(i)*dt
END DO
amsum = 0.0
i = 1
DO WHILE (amsum < 0.5*fwhmsum)
  amsum = amsum + fadeco(i)*dt
  fwhmtime = tt(i)
  i = i + 1
END DO
!! WRITE(*,*) 'FWHM range and integral:', fwhm1,fwhm2,fwhmsum,fwhmtime
  
sisum = SUM(fadeco)*dt

! Background level
bgsumde = 0.0
nbgde = 0
DO i = 1, nt
!  IF (tt(i) > 10.0 .AND. tt (i) < 90.0) THEN
!    bgsumde = bgsumde + fadeco(i)
!    nbgde = nbgde + 1
!  END IF
  IF (tt(i) > 450.0 .AND. tt (i) < 950.0) THEN
    bgsumde = bgsumde + fadeco(i)
    nbgde = nbgde + 1
  END IF
END DO

bgsumde = bgsumde/REAL(nbgde)

WRITE(19,*) 'Background level (deconvolved):', bgsumde
!! WRITE(*,*) 'Background level (deconvolved):', bgsumde

! Reconstructed number of photoelectrons
repe = (sampli-bgsumde)/refeped
repede = fwhmsum/refepedsum

!! WRITE(*,*) 'Deconvolved signal, amplitude:', sampli, sisum, sampsum, resum
!! WRITE(*,*) 'Deconvolved signal, arrival time:', recoti,ptime(1)
!! WRITE(*,*) 'Deconvolved signal, number of photoelectrons:', repe,repede

WRITE(19,*) 'Deconvolved signal, FWHM range and integral:', fwhm1,fwhm2,fwhmsum,fwhmtime
WRITE(19,*) 'Deconvolved signal, amplitude:', sampli, sisum, sampsum, resum
WRITE(19,*) 'Deconvolved signal, arrival time:', recoti,ptime(1)
WRITE(19,*) 'Deconvolved signal, number of photoelectrons:', repe,repede
WRITE(19,*) 'Deconvolved signal, fine amplitude and arrival time:', samplide,decotin
!WRITE(19,*) 'FADC delay time:', tadc

! Shaped signal: three moving averages of the deconvolved signal

! 250 MHz
tave = 2.0*tsam
! 1 GHz
!tave = 8.0*tsam
!tave = 8.0

CALL herin(nt,ni,tt,ti,fadeco,fadin)
CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

CALL herin(nt,ni,tt,ti,fasa,fadin)
CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

CALL herin(nt,ni,tt,ti,fasa,fadin)
CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

CALL herin(nt,ni,tt,ti,fasa,fadin)
CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

! Two more for 1 GHz
!CALL herin(nt,ni,tt,ti,fasa,fadin)
!CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

!CALL herin(nt,ni,tt,ti,fasa,fadin)
!CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

sampli = -100000.0
DO i = 1, nt
  IF (tt(i) > 420.0) EXIT
  IF (fasa(i) > sampli) THEN
    sampli = fasa(i)
    recoti = tt(i)
  END IF
END DO

sisum = SUM(fasa)*dt

DO i = 1, nt
! fadc.shaped
  WRITE(35,*) tt(i),fasa(i)
END DO
CLOSE(35)


! Background level
bgsumsh = 0.0
nbgsh = 0
DO i = 1, nt
!  IF (tt(i) > 10.0 .AND. tt (i) < 80.0) THEN
!    bgsumsh = bgsumsh + fasa(i)
!    nbgsh = nbgsh + 1
!  END IF
  IF (tt(i) > 450.0 .AND. tt (i) < 950.0) THEN
    bgsumsh = bgsumsh + fasa(i)
    nbgsh = nbgsh + 1
  END IF
END DO

bgsumsh = bgsumsh/REAL(nbgsh)

WRITE(19,*) 'Background level (shaped):', bgsumsh
!! WRITE(*,*) 'Background level (shaped):', bgsumsh

! Fine interpolation for a better time resolution
!CALL linin(nt,ni,tt,ti,fasa,fasin)
CALL herin(nt,ni,tt,ti,fasa,fasin)

samplin = -100000.0
DO i = 1, ni
  IF (ti(i) > 420.0) EXIT
  IF (fasin(i) > samplin) THEN
    samplin = fasin(i)
    recotin = ti(i)
  END IF
END DO

DO i = 1, ni
! fadc.shapin
  WRITE(36,*) ti(i),fasin(i)
END DO
CLOSE(36)

! Reconstructed number of photoelectrons
repe = (sampli-bgsumsh)/refepes

!! WRITE(*,*) 'Shaped signal, amplitude:', sampli, sisum, sampsum, resum
!! WRITE(*,*) 'Shaped signal, arrival time:', recoti,ptime(1)
!! WRITE(*,*) 'Shaped signal, number of photoelectrons:', repe
!! WRITE(*,*) 'Shaped signal, fine amplitude and arrival time:', samplin,recotin

WRITE(19,*) 'Shaped signal, amplitude:', sampli, sisum, sampsum, resum
WRITE(19,*) 'Shaped signal, arrival time:', recoti,ptime(1)
WRITE(19,*) 'Shaped signal, number of photoelectrons:', repe
WRITE(19,*) 'Shaped signal, fine amplitude and arrival time:', samplin,recotin
WRITE(19,*) 'FADC delay time:', tadc


CONTAINS


! NSB arrival times
SUBROUTINE nsbat(nsbra,tsta,tsto,nsbnu,nsbti)

REAL, INTENT(IN) :: nsbra,tsta,tsto
INTEGER, INTENT(OUT) :: nsbnu
REAL, DIMENSION(100000), INTENT(OUT) :: nsbti

REAL :: raa,rab,lambda2,rate
REAL(KIND=8) :: lambda,nsbnur,wran,pran
INTEGER :: i

! Find Poisson distributed number of NSB arrivals
! NSB rate converted to photons per ns
rate = (tsto-tsta)*nsbra/1000.0
lambda = REAL(rate,KIND=8)
lambda2 = 6.0*rate
DO
  CALL RANDOM_NUMBER(raa)
  nsbnu = NINT(lambda2*raa)
  nsbnur = REAL(nsbnu,KIND=8)
  wran = dpois(nsbnur,lambda)
  pran = DEXP(wran)
  CALL RANDOM_NUMBER(rab)
  IF (rab < pran) EXIT
END DO

! Simulate arrival times
DO i = 1, nsbnu
  CALL RANDOM_NUMBER(raa)
  nsbti(i) = tsta + raa*(tsto-tsta)
END DO

RETURN

END SUBROUTINE nsbat

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

! Linear interpolation
SUBROUTINE linin(ni,no,xi,xo,fi,fo)

REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL :: dxi,dxo,ai,bi

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3

dxi = 0.5*(xi(2)-xi(1))
dxo = 0.5*(xo(2)-xo(1))

DO i = 1, no
  fo(i) = 0.0
! Select new bins covering input bins
  IF ( xo(i)-dxo .GE. xi(1)+dxi .AND. xo(i)-dxo .LE. xi(ni)+dxi ) THEN
    DO j = 1, ni
      IF ( xi(j) <= xo(i) .AND. xi(j+1) > xo(i) ) THEN
        j1 = j
        j2 = j+1
      END IF
    END DO
    ai = (fi(j2)-fi(j1))/(xi(j2)-xi(j1))
    bi = (xi(j2)*fi(j1)-xi(j1)*fi(j2))/(xi(j2)-xi(j1))
    fo(i) = ai*xo(i) + bi
  END IF
END DO
RETURN  
END SUBROUTINE linin
      
! Parabolic interpolation
SUBROUTINE parin(ni,no,xi,xo,fi,fo)

REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL :: dxi,dxo,ai,bi

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3

dxi = 0.5*(xi(2)-xi(1))
dxo = 0.5*(xo(2)-xo(1))

DO i = 1, no
  fo(i) = 0.0
! Select new bins covering input bins
  IF ( xo(i)-dxo .GE. xi(1)+dxi .AND. xo(i)-dxo .LE. xi(ni)+dxi ) THEN
    DO j = 1, ni
      det1 = ABS(xi(j)-xo(i))
      det2 = ABS(xi(j+1)-xo(i))
      IF (xi(j) <= xo(i) .AND. xi(j+1) > xo(i)) THEN
        j1 = j
        j2 = j+1
        j3 = j+2
        IF (det1 < det2) THEN
          j1 = j-1
          j2 = j
          j3 = j+1
        END IF
      END IF
    END DO
    fo(i) = parabola(xo(i),xi(j1),xi(j2),xi(j3),fi(j1),fi(j2),fi(j3))
  END IF
END DO
RETURN  
END SUBROUTINE parin

! Hermite interpolation
SUBROUTINE herin(ni,no,xi,xo,fi,fo)

REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL(KIND=8) :: bias,tension,x1,x2,x3,x4,y1,y2,y3,y4,dxi,dxo,xo2,xo3,yy1,yy2

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,j4

tension = 0.0
bias = 0.0

dxi = 0.5*(xi(2)-xi(1))
dxo = 0.5*(xo(2)-xo(1))

DO i = 1, no
  fo(i) = 0.0
  IF ( xo(i) .GE. xi(2) .AND. xo(i) .LE. xi(ni-1) ) THEN
    DO j = 2, ni-1
      det1 = ABS(xi(j)-xo(i))
      det2 = ABS(xi(j+1)-xo(i))
      IF (xi(j) <= xo(i) .AND. xi(j+1) > xo(i)) THEN
        j1 = j-1
        j2 = j
        j3 = j+1
        j4 = j+2
      END IF
    END DO
    xo2 = (xo(i)-xi(j2))*(xo(i)-xi(j2))/((xi(j3)-xi(j2))*(xi(j3)-xi(j2)))
    xo3 = (xo(i)-xi(j2))*(xo(i)-xi(j2))*(xo(i)-xi(j2))/((xi(j3)-xi(j2))*(xi(j3)-xi(j2))*(xi(j3)-xi(j2)))
    yy1 = (fi(j2)-fi(j1))*(1.0+bias)*(1.0-tension)/2.0
    yy1 = yy1 + (fi(j3)-fi(j2))*(1.0-bias)*(1.0-tension)/2.0
    yy2 = (fi(j3)-fi(j2))*(1.0+bias)*(1.0-tension)/2.0
    yy2 = yy2 + (fi(j4)-fi(j3))*(1.0-bias)*(1.0-tension)/2.0
    x1 = 2.0*xo3 - 3.0*xo2 + 1.0
    x2 = xo3 - 2.0*xo2 + (xo(i)-xi(j2))/(xi(j3)-xi(j2))
    x3 = xo3 - xo2
    x4 = -2.0*xo3 + 3.0*xo2
    fo(i) = x1*fi(j2) + x2*yy1 + x3*yy2 + x4*fi(j3)
    IF (xo(i) > 105.0 .AND. xo(i) < 106.0) THEN
!      WRITE(*,*) 'herin',xo(i),yy1,yy2,x1,x2,x3,x4,fo(i)
    END IF
  END IF
END DO
RETURN
END SUBROUTINE herin

! /*
!   Tension: 1 is high, 0 normal, -1 is low
!   Bias: 0 is even,
!         positive is towards first segment,
!         negative towards the other
!*/
!double HermiteInterpolate(
!   double y0,double y1,
!   double y2,double y3,
!   double mu,
!   double tension,
!   double bias)
!{
!   double m0,m1,mu2,mu3;
!   double a0,a1,a2,a3;

!	mu2 = mu * mu;
!	mu3 = mu2 * mu;
!   m0  = (y1-y0)*(1+bias)*(1-tension)/2;
!   m0 += (y2-y1)*(1-bias)*(1-tension)/2;
!   m1  = (y2-y1)*(1+bias)*(1-tension)/2;
!   m1 += (y3-y2)*(1-bias)*(1-tension)/2;
!   a0 =  2*mu3 - 3*mu2 + 1;
!   a1 =    mu3 - 2*mu2 + mu;
!   a2 =    mu3 -   mu2;
!   a3 = -2*mu3 + 3*mu2;

!   return(a0*y1+a1*m0+a2*m1+a3*y2);
!}

! Cubic spline interpolation
SUBROUTINE curin(ni,no,xi,xo,fi,fo)
REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL :: bias,tension,x1,x2,x3,x4,y1,y2,y3,y4,dxi,dxo,xo2,xo3,yy1,yy2

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,j4

tension = 0.0
bias = 0.0

dxi = 0.5*(xi(2)-xi(1))
dxo = 0.5*(xo(2)-xo(1))

DO i = 1, no
  fo(i) = 0.0
  IF ( xo(i) .GE. xi(2) .AND. xo(i) .LE. xi(ni-1) ) THEN
    DO j = 2, ni-1
      det1 = ABS(xi(j)-xo(i))
      det2 = ABS(xi(j+1)-xo(i))
      IF (xi(j) <= xo(i) .AND. xi(j+1) > xo(i)) THEN
        j1 = j-1
        j2 = j
        j3 = j+1
        j4 = j+2
      END IF
    END DO
    xo2 = (xo(i)-xi(j2))*(xo(i)-xi(j2))/((xi(j3)-xi(j2))*(xi(j3)-xi(j2)))
    xo3 = (xo(i)-xi(j2))*(xo(i)-xi(j2))*(xo(i)-xi(j2))/((xi(j3)-xi(j2))*(xi(j3)-xi(j2))*(xi(j3)-xi(j2)))
    x1 = -0.5*fi(j1) + 1.5*fi(j2) - 1.5*fi(j3) + 0.5*fi(j4)
    x2 = fi(j1) - 2.5*fi(j2) + 2.0*fi(j3) - 0.5*fi(j4)
    x3 = -0.5*fi(j1) + 0.5*fi(j3)
    x4 = fi(j2)
    fo(i) = x1*xo3 + x2*xo2 + x3*(xo(i)-xi(j2))/(xi(j3)-xi(j2)) + x4
    IF (xo(i) > 105.0 .AND. xo(i) < 106.0) THEN
!      WRITE(*,*) 'curin',xo(i),x1,x2,x3,x4,fo(i)
    END IF
  END IF
END DO
RETURN  
END SUBROUTINE curin

!double CubicInterpolate(
!   double y0,double y1,
!   double y2,double y3,
!   double mu)
!{
!   double a0,a1,a2,a3,mu2;

!   mu2 = mu*mu;
!   a0 = y3 - y2 - y0 + y1;
!   a1 = y0 - y1 - a0;
!   a2 = y2 - y0;
!   a3 = y1;

!   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
!}

!Paul Breeuwsma proposes the following coefficients for a smoother interpolated curve, 
!which uses the slope between the previous point and the next as the derivative at the 
!current point. This results in what are generally referred to as Catmull-Rom splines.

!   a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
!   a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3;
!   a2 = -0.5*y0 + 0.5*y2;
!   a3 = y1;

! Polynomial interpolation, Numerical Recipes
SUBROUTINE porin(ni,no,xi,xo,fi,fo)
REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL, DIMENSION(10) :: xp,yp
REAL :: bias,tension,x1,x2,x3,x4,y1,y2,y3,y4,dxi,dxo,xo2,xo3,yy1,yy2,dyp

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,j4,k

tension = 0.0
bias = 0.0

dxi = 0.5*(xi(2)-xi(1))
dxo = 0.5*(xo(2)-xo(1))

DO i = 1, no
  fo(i) = 0.0
  IF ( xo(i) .GE. xi(2) .AND. xo(i) .LE. xi(ni-1) ) THEN
    DO j = 2, ni-1
      det1 = ABS(xi(j)-xo(i))
      det2 = ABS(xi(j+1)-xo(i))
      IF (xi(j) <= xo(i) .AND. xi(j+1) > xo(i)) THEN
        j1 = j-1
        j2 = j
        j3 = j+1
        j4 = j+2
      END IF
    END DO
    xp(1) = xi(j1)
    xp(2) = xi(j2)
    xp(3) = xi(j3)
    xp(4) = xi(j4)
    yp(1) = fi(j1)
    yp(2) = fi(j2)
    yp(3) = fi(j3)
    yp(4) = fi(j4)
    CALL polint(xp,yp,4,xo(i),fo(i),dyp)
  END IF
END DO
RETURN  
END SUBROUTINE porin

SUBROUTINE polint(xa,ya,n,x,y,dy)
INTEGER :: n,nmax
REAL :: dy,x,y,xa(n),ya(n)
! Largest anticipated value of n
PARAMETER (nmax=10)
! Given arrays xa and ya, each of length n, and given a value x, this routine
! returns value y, and an error estimate dy. If P(x) is the polynomial of degree
! N-1 such that P(xa)=ya, i=1,...,n then the returned value y = P(x)
INTEGER :: i,m,ns
REAL :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)
ns = 1
dif = ABS(x-xa(1))
DO i=1,n
  dift = ABS(x-xa(i))
  IF (dift .LT. dif) THEN
    ns = i
    dif = dift
  END IF
  c(i) = ya(i)
  d(i) = ya(i)
END DO
y = ya(ns)
ns = ns-1
DO m=1,n-1
  DO i=1,n-m
    ho = xa(i)-x
    hp = xa(i+m)-x
    w = c(i+1)-d(i)
    den = ho - hp
    IF (den .EQ. 0.0) WRITE(*,*) 'failure in polint', x,xa,m,i,ho,hp
    den = w/den
    d(i) = hp*den
    c(i) = ho*den
  END DO
  IF (2*ns .LT. n-m) THEN
    dy = c(ns+1)
  ELSE
    dy = d(ns)
    ns = ns-1
  END IF
  y = y + dy
END DO
RETURN
END SUBROUTINE polint  
     
! Moving average 
SUBROUTINE movave(ni,no,xi,xo,fi,fo,ts)
REAL, DIMENSION(10000), INTENT(IN) :: xi,fi,xo
REAL, DIMENSION(10000), INTENT(OUT) :: fo
REAL, INTENT(IN) :: ts
REAL :: asum,wsum,dfo,dxi,we1,we2,we3,we4

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,dj

dxi = (xi(ni)-xi(1))/REAL(ni-1)

DO i = 1, no
  fo(i) = 0.0
  IF (xo(i)-0.5*ts > xi(1) .AND. xo(i)+0.5*ts < xi(ni)) THEN
    DO j = 1, ni
      IF (xo(i)-0.5*ts > xi(j) .AND. xo(i)-0.5*ts <= xi(j+1)) j1 = j+1
      IF (xo(i)+0.5*ts >= xi(j) .AND. xo(i)+0.5*ts < xi(j+1)) j2 = j
    END DO
    asum = 0.0
    wsum = 0.0
    we1 = 0.0
    we2 = 0.0
    we3 = 0.0
    we4 = 0.0
! Partial overlap 
    IF (xo(i)-0.5*ts < xi(j1-1)+0.5*dxi) THEN
      we1 = (xi(j1-1)+0.5*dxi-xo(i)+0.5*ts)/dxi
      asum = asum + we1*fi(j1-1) + fi(j1)
      wsum = wsum + we1 + 1.0
    END IF
    IF (xo(i)-0.5*ts >= xi(j1-1)+0.5*dxi) THEN
      we2 = (xi(j1)+0.5*dxi-xo(i)+0.5*ts)/dxi
      asum = asum + we2*fi(j1)
      wsum = wsum + we2
    END IF
    DO j = j1+1, j2-1
      asum = asum + fi(j)
      wsum = wsum + 1.0
    END DO
    IF (xo(i)+0.5*ts > xi(j2+1)-0.5*dxi) THEN
      we3 = (xi(j2+1)+0.5*dxi-xo(i)-0.5*ts)/dxi
      asum = asum + we3*fi(j2+1) + fi(j2)
      wsum = wsum + we3 + 1.0
    END IF
    IF (xo(i)+0.5*ts <= xi(j2)+0.5*dxi) THEN
      we4 = (xi(j2)+0.5*dxi-xo(i)-0.5*ts)/dxi
      asum = asum + we4*fi(j2)
      wsum = wsum + we4
    END IF
    dj = j2-j1+1
    fo(i) = asum/wsum
!    fo(i) = asum/REAL(j2-j1+1)
!    dfo = fo(i) - fo(i-1)
!    IF (xo(i) > 102.4 .AND. xo(i) < 103.6) WRITE(*,*) 'movave:', dxi,xo(i),fo(i),dj,dfo,we1,we2,we3,we4,wsum
  END IF
END DO
RETURN
END SUBROUTINE movave

! Hybrid gaussian
FUNCTION hyga(x,am,t0,si,tu)
  REAL, INTENT(IN) :: x,am,t0,si,tu
  REAL :: hyga,sih
  sih = 2.0*si*si + (x-t0)*tu
  hyga = am*EXP(-(x-t0)*(x-t0)/sih)
  RETURN
END FUNCTION hyga

! Two exponents
FUNCTION dexe(x,am,t0,tar,tad)
  REAL, INTENT(IN) :: x,am,t0,tar,tad
  REAL :: dexe
  dexe = am*(1.0-EXP(-(x-t0)/tar))*EXP(-(x-t0)/tad)
  RETURN
END FUNCTION dexe

! Random number for a uniform (box) distribution over some range
FUNCTION bodi(bce,bra)
  REAL, INTENT(IN) :: bce,bra
  REAL :: bodi,brn
  CALL RANDOM_NUMBER(brn)
  bodi = bce + (brn-0.5)*bra
  RETURN
END FUNCTION bodi

! Random amplitude for a given SPE distribution
FUNCTION speresp(as,spx,spr,ns)
  REAL, DIMENSION(10000), INTENT(IN) :: spx,spr
  REAL, INTENT(IN) :: as
  INTEGER, INTENT(IN) :: ns
  REAL :: speresp,ra,ram,rasp
  INTEGER :: ira
  DO
    CALL RANDOM_NUMBER(ra)
    ram = spx(1) + ra*(spx(ns)-spx(1))
    ira = ANINT((ram-spx(1))/(spx(2)-spx(1)))
    rasp = spr(ira)
    CALL RANDOM_NUMBER(ra)
    IF (ra < rasp) EXIT
  END DO
  speresp = ram*as
  RETURN
END FUNCTION speresp

! Random number for a normal distribution
FUNCTION rano(noc,nos)
  REAL, INTENT(IN) :: noc,nos
  REAL :: rano,nra,dino
  DO
    CALL RANDOM_NUMBER(nra)
    rano = noc + 12.0*(nra-0.5)*nos
    dino = EXP(-(rano-noc)*(rano-noc)/(2.0*nos*nos))
    CALL RANDOM_NUMBER(nra)
    IF (nra < dino) EXIT
  END DO
  RETURN
END FUNCTION rano

! Convolution with the amplifier/filters response
FUNCTION amfire(ipu,pud,tid,ama,amb,amt)
  REAL, INTENT(IN) :: ama,amb,amt
  REAL, DIMENSION(10000), INTENT(IN) :: pud,tid
  REAL :: amfire,ares,sumco,detu,dett
  INTEGER, INTENT(IN) :: ipu
  INTEGER :: i
  sumco = 0.0
  DO i = 1, ipu
    detu = tt(i)-tt(1)
    ares = (ama+amb*detu)*exp(-detu/amt)
    sumco = sumco + pud(ipu-i)*ares*(tt(2)-tt(1))
  END DO
  amfire = sumco
  RETURN
END FUNCTION amfire

! Cubic spline
FUNCTION parabola(x,y1,y2,y3,f1,f2,f3)   
  REAL, INTENT(IN) :: x,y1,y2,y3,f1,f2,f3
  REAL :: parabola,a,b,c
  REAL :: dy1,dy2,dyy1,dyy2,df1,df2,atest
  dy1 = y2 - y1
  dy2 = y3 - y1
  dyy1 = y2*y2 - y1*y1
  dyy2 = y1*y1 - y3*y3
  df1 = f2 - f1
  df2 = f1 - f3
!  IF (x > 99.0 .AND. x < 100.0) WRITE(*,*) 'parabola', dy1,dy2,dyy1,dyy2,df1,df2
  
  a = (dy2*df1+dy1*df2)/(dyy1*dy2+dyy2*dy1)
!  IF (ISNAN(a)) THEN
!    atest = dyy1*dy2+dyy2*dy1
!    WRITE(*,*) 'dy1,dy2,dyy1,dyy2', dy1,dy2,dyy1,dyy2,atest
!  END IF
  b = (df1-a*dyy1)/(dy1)    
  c = f1 - a*y1*y1 - b*y1  
!  dy1 = a*(y2+0.1)*(y2+0.1) + b*(y2+0.1) + c 
  parabola = a*x*x + b*x + c
!  IF (ISNAN(parabola)) THEN
!    atest = dyy1*dy2+dyy2*dy1
!    WRITE(*,*) 'dy1,dy2,dyy1,dyy2', dy1,dy2,dyy1,dyy2,atest
!  END IF
!  IF (x > 70.0 .AND. x < 71.0) WRITE(*,*) 'parabola', x,y2,a,b,c,parabola,dy1
  RETURN
END FUNCTION parabola
  
FUNCTION dpois(x,lb)
      
  REAL(KIND=8), INTENT(IN) :: lb,x
  REAL(KIND=8) :: dpois
        
  REAL(KIND=8) :: pi2
      
  pi2 = 6.283185307179586476925286

  IF (lb.EQ.0.0D0) THEN
      
    IF (x.EQ.0.0D0) dpois = 1.0D0
    IF (x.NE.0.0D0) dpois = 0.0D0

  RETURN
	
  END IF
      
  IF (x.EQ.0.0D0) THEN
      
!    dpois = EXP(-lb)
    dpois = -lb
   
    RETURN
	
  END IF
       
  IF (x.LT.16.0) THEN 
    
!    dpois = EXP(-stor(x)-bd0(x,lb))/SQRT(pi2*x)
    dpois = -stor(x)-bd0(x,lb) - DLOG(DSQRT(pi2*x))

  END IF	
      
  IF (x.GE.16.0) THEN 
   
!    dpois = EXP(-stir(x)-bd0(x,lb))/SQRT(pi2*x)
    dpois = -stir(x)-bd0(x,lb) - DLOG(DSQRT(pi2*x))

  END IF	
      
  RETURN
  
  END FUNCTION dpois       	
      
      
  FUNCTION stir(z)
      
  REAL(KIND=8), INTENT(IN) :: z
  REAL(KIND=8) :: stir
      
  REAL(KIND=8) :: s0, s1, s2, s3, s4, zz, z2
  REAL(KIND=8), DIMENSION(16) :: sfe(16)
  INTEGER :: iz
      
      s0 = 0.08333333333333333333333D0
      s1 = 0.0027777777777777777777778D0
      s2 = 0.00079365079365079365079365D0
      s3 = 0.000595238095238095238095238D0
      s4 = 0.0008417508417508417508417508D0
      
      sfe(1) = 0.0D0
      sfe(2) = 0.081061466795327258219670264
      sfe(3) = 0.041340695955409294093822081
      sfe(4) = 0.0276779256849983391487892927
      sfe(5) = 0.020790672103765093111522771
      sfe(6) = 0.0166446911898211921631948653
      sfe(7) = 0.013876128823070747998745727
      sfe(8) = 0.0118967099458917700950557241
      sfe(9) = 0.010411265261972096497478567
      sfe(10) = 0.0092554621827127329177286366
      sfe(11) = 0.008330563433362871256469318
      sfe(12) = 0.0075736754879518407949720242
      sfe(13) = 0.006942840107209529865664152
      sfe(14) = 0.0064089941880042070684396310
      sfe(15) = 0.005951370112758847735624416
      sfe(16) = 0.0055547335519628013710386899   
 
      iz = INT(z)
 
      z2 = z*z
      
      IF (iz.LT.16) THEN
      
        stir = sfe(iz+1)
	
        RETURN
          	
      END IF
      
      IF (iz.GT.500) THEN
      
        stir = (s0-s1/z2)/z
	
	RETURN
	
      END IF
      
      IF (iz.GT.80) THEN
      
        stir = (s0-(s1-s2/z2)/z2)/z
	
	RETURN
	
  END IF
      
  IF (iz.GT.35) THEN
      
    stir = (s0-(s1-(s2-s3/z2)/z2)/z2)/z
	
    RETURN
	
  END IF
      
  stir = (s0-(s1-(s2-(s3-s4/z2)/z2)/z2)/z2)/z
      
  RETURN
      
END FUNCTION stir
      
      
FUNCTION bd0(x,np)
      
REAL(KIND=8), INTENT(IN) :: x,np
REAL(KIND=8) :: bd0
  
REAL(KIND=8) :: ej, s, s1, v, dx, d1, d2

INTEGER j
      
      d1 = REAL(DABS(x-np),KIND=8)
      d2 = 0.1d0*REAL(DABS(x+np),KIND=8)
      
      IF (d1.LT.d2) THEN
      
        s = (x-np)*(x-np)/(x+np)
        v = (x-np)/(x+np)
        ej = 2.0D0*x*v
	v = v*v
        DO j = 1, 1000
          ej = ej*v
	  s1 = s+ej/(2.0D0*j+1.0D0)
          IF (s1.EQ.s) THEN
	    bd0 = s1
	    RETURN
	  END IF
          s = s1
	END DO 
      END IF
      bd0 = x*REAL(DLOG(x/np),KIND=8)+np-x
      RETURN
      
END FUNCTION bd0      
      
      
FUNCTION stor(z)
      
REAL(KIND=8), INTENT(IN) :: z
REAL(KIND=8) :: stor

REAL(KIND=8) :: pi2, gz, as
      
      pi2 = 6.283185307179586476925286D0

      gz = DEXP(gammln(z))
      
      as = DSQRT(pi2*z)
      
      stor = DLOG(gz*DEXP(z)/(as*z**(z-1.0D0)))
      
  RETURN
  
END FUNCTION stor
      

FUNCTION gamdis(el,ek,en)

REAL(KIND=8), INTENT(IN) :: el,ek,en
REAL(KIND=8) :: gamdis
REAL(KIND=8) :: enor
      
IF (ek.LE.0.0D0) THEN
  gamdis = DEXP(-el*en)*en**(-1.0D0)
  RETURN
END IF	
      
enor = el**ek/DEXP(gammln(ek))

!WRITE(*,*) 'enor', enor
      
gamdis = enor*DEXP(-el*en)*en**(ek-1.0D0)
      
RETURN

END FUNCTION gamdis
      
        	
      
FUNCTION gammln(xx) 
!     Returns the value ln[\Gamma (xx)] for xx > 0
      
REAL(KIND=8), INTENT(IN) :: xx
REAL(KIND=8) :: gammln

INTEGER :: j 
REAL(KIND=8) :: ser,stp,tmp,x,y
REAL(KIND=8), DIMENSION(6) :: cof

  cof(1) = 76.18009172947146D0
  cof(2) = -86.50532032941677D0
  cof(3) = 24.01409824083091D0
  cof(4) = -1.231739572450155D0
  cof(5) = .1208650973866179D-2
  cof(6) = -.5395239384953D-5
  stp = 2.5066282746310005D0

x=xx 
y=x 
tmp=x+5.5D0 
tmp=(x+0.5D0)*DLOG(tmp)-tmp 
ser=1.000000000190015D0 

DO j = 1, 6
  y = y + 1.0D0 
  ser = ser + cof(j)/y 
END DO 
      
gammln=tmp+DLOG(stp*ser/x)

RETURN 

END FUNCTION gammln


END PROGRAM flasim
