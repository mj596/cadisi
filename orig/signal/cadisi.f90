PROGRAM cadisi

! Version of 31 July 2013 (cleaned and faster version of flasim)
! Simulates/computes digital camera (FlashCam) signal reconstruction
! Both photon detector options: GAPD and photomultiplier
! Option to use pedestal computed from the previous (calibration) run
! instead the internal pedestal (signal level outside the Cherenkov pulse)


IMPLICIT NONE

REAL, DIMENSION(200000) :: t1,t2,tt,te,pude,hit,pufi,ta,fadi,ti,fadan,fadepi
REAL, DIMENSION(200000) :: fada,tu,fadu,fade,fades,fadp,fadep,fadesa,fadar,fadeco
REAL, DIMENSION(200000) :: fady,fasa,fast,fadc,fadin,wnoise,fadn,fasin,fadeni,fadb
REAL, DIMENSION(200000) :: nsbart,cerart,atime,ptime,spex,sper,am,fader,fadap,fadapi
REAL, DIMENSION(200000) :: tenart,tfine,afine,pudefine,pufine
REAL :: atimemin,atimemax,sumati,artime,detime,sitime,time0,uh,th,sh,ah,ytest
REAL :: desum,amsum,ares,deta,tsam,tstart,tadc,dt,onepe,tad,fasum,ai,tic,deti
REAL :: awe,starti,sisum,det1,det2,bi,detimin,tam1,tam2,dti,swe,resum,sampli
REAL :: recoti,sampsum,repe,refepe,pudemax,pufisum,pufimax,sampmax,isamsum,fisum
REAL :: isammax,misasum,misamax,tave,dth,af1,af2,af3,wnom,wnos,wnomean,wnoisa
REAL :: tno,nsbrate,tstop,tra,tde,fwhm1,fwhm2,fwhmsum,fwhmtime,refeped,refepes
REAL :: refepedsum,repede,bgsumde,bgsumsh,spet,samplin,recotin,recota,recotu
REAL :: decotin,samplide,posisum,nesisum,pomean,nemean,anesisum,anemean
REAL :: dif1,dif2,dif3,dif4,fasumi,defa,signorm,onepa,wnomeanpe,tenrate
REAL :: deltime,tjitter,tsam14,tsam34,tsam12,dt12,bi1,bi2,bo1,bo2,punor
REAL :: fasamax,tisamax,maxamp,maxtim,racer,ramean,rasigma,ravar
REAL :: phadc,dtfine,tamp,wfac,tfines,pedelevel
REAL :: fanoise,pevar,pepar1,pepar2,pepar3,qevar,gavar,sevar

INTEGER, DIMENSION(200000) :: nf,nno,pocer,ponu
INTEGER, DIMENSION(2000) :: me
INTEGER, DIMENSION(1) :: imax
INTEGER :: i,k,j,nc,istart,istop,ne,nseed,j0,j1,j2,j3,jj1,jj2,jj3,na,ii,ni,nt
INTEGER :: nsbnum,np,ihm1,ihm2,nbgde,nbgsh,nspe,npos,nneg,ncer,tennum,ncerpo,ide
INTEGER :: pomin,pomax,nfine,ifine0,ifine1,ifine2,ipor,ispe,itim,ipr,ipha,ifines
INTEGER, DIMENSION(:), ALLOCATABLE :: seed

REAL :: cerate,nsrate,ton,toff,events,photons,pi,sigma1p,gno,gen,gex,celec1p
REAL :: sigma0,celec0,gex0,gno0,gen0,darate,sigma1,celec1,gnop,gexp,genp,sigmap
REAL :: alpha,we,ala,beta,era,alpa,celecn,sigman,ar1,ar2,sia,tevents,wep
REAL :: sump,sumpe,mumax,pmax,mumean,cemu,sumap,mumedian,sigmamu,clower,cupper
REAL :: sigmal,sigmau,sumpi,sumpa,mimean,sigmami,mimedian,nsmean,nssigma
REAL :: cnsb,wnsb,npnsb,sumns,sumfadis,ranu,rana,maxfadis,den0,dem0
REAL :: sigma02,sigma12,telaps1,telaps2
REAL, DIMENSION(2000) :: nfoto,pnsb
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

! Read input parameters
OPEN(1,FILE='cadisi.par')
! Cherenkov photons
READ(1,*) ncer
! Randomization flag for the mean number of Cherenkov photons
READ(1,*) ipor
! Mean arrival time of Cherenkov photons
READ(1,*) time0
! Flag for randomization of arrival time of Cherenkov photons
READ(1,*) itim
! Time spread of Cherenkov photons, box distribution (ns)
READ(1,*) deltime
! Time jitter for photon detector
READ(1,*) tjitter
! Random phase flag
READ(1,*) ipha
! ADC sampling period
READ(1,*) tsam
! Flag for SPE 
READ(1,*) ispe
! GAPD/PMT amplitude corresponding to 1 photoelectron
READ(1,*) onepe
! Amplifier amplitude corresponding to 1 photoelectron, used only for electronic noise
READ(1,*) onepa
! White noise fractio factor
READ(1,*) wfac
! White noise mean amplitude
READ(1,*) wnom
! NSB rate MHz
READ(1,*) nsbrate
! Thermal noise of GAPD, 6 MHz
READ(1,*) tenrate
! Derivative added to the averaged signal with factor defa
READ(1,*) defa
! Single photon signal normalization
READ(1,*) signorm
! Reference amplitude of deconvolved signal
READ(1,*) refeped
! Reference amplitude of shaped signal
READ(1,*) refepes
! Start time of computed period (ns)
READ(1,*) tstart
! Stop time of computed period (ns)
READ(1,*) tstop
! Coarse time bin (ns)
READ(1,*) dt
! Fine time bin (ns)
READ(1,*) dti
! Finest time bin (ns)
READ(1,*) dtfine
! Computation limits for background level (pedestal)
READ(1,*) bi1
READ(1,*) bo1
READ(1,*) bi2
READ(1,*) bo2
! Printing flag
READ(1,*) ipr
! Pedstal level ( > 0 - use calibration run instead computing sumbg)
READ(1,*) pedelevel
! FADC noise
READ(1,*) fanoise
! pedestal variations (e.g. temperature)
READ(1,*) pevar
! channel to channel pedestal varations
READ(1,*) pepar1
! FADC deviations for a single channel
READ(1,*) pepar2 
! error of initial calibration
READ(1,*) pepar3
! quantum efficiency variations
READ(1,*) qevar
! gain variations
READ(1,*) gavar
! sensitivity variations
READ(1,*) sevar

CLOSE(1)

nt = NINT((tstop-tstart)/dt)
nfine = NINT((tstop-tstart)/dtfine)
IF (ipr == 1) WRITE(*,*) 'Coarse and finest interpolation bins nt nfine:', nt, nfine

! White noise standard deviation, Cracow: 0.15 mV noise for 1 mV of one p.e. signal
wnos = wfac*onepa

tsam12 = 0.5*tsam
tsam14 = 0.25*tsam
tsam34 = 0.75*tsam

nseed = 17
CALL RANDOM_SEED(size = nseed)
ALLOCATE(seed(nseed))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i, i = 1, nseed) /)
CALL RANDOM_SEED(PUT = seed)
DEALLOCATE(seed)

!! WRITE(*,*) 'Coarse and fine bins:', nt, ni

! Single photoelectron spectrum with cross-talk and afterpulses
! fitted to the mppc_uknown.dat data

OPEN(10,FILE='cadisi.results')
OPEN(11,FILE='mppcspenormean.dat')
!OPEN(12,FILE='speresp.test')

IF (ipr == 1) THEN
  OPEN(21,FILE='cherenkov.photons')
  OPEN(22,FILE='nsb.photons')
  OPEN(23,FILE='thermal.photons')
  OPEN(24,FILE='detector.signals')
  OPEN(31,FILE='detector.pulse')
  OPEN(32,FILE='amplifier.pulse')
  OPEN(33,FILE='fadc.samples')
!  OPEN(34,FILE='fadc.inter')
  OPEN(35,FILE='fadc.aver')
!  OPEN(36,FILE='fadc.averin')
  OPEN(37,FILE='fadc.dave')
!  OPEN(38,FILE='fadc.daves')
  OPEN(39,FILE='fadc.daver')
!  OPEN(40,FILE='fadc.davesa')
  OPEN(41,FILE='fadc.deco')
  OPEN(42,FILE='fadc.shaped')
  OPEN(43,FILE='fadc.shapin')
!  OPEN(44,FILE='fadc.sampled')
  OPEN(51,FILE='flasim.test')
  OPEN(52,FILE='flasim.dist')
  OPEN(53,FILE='photon.bunches')
!  OPEN(54,FILE='white.noise')
!  OPEN(55,FILE='random.photons')
END IF

nspe = 0
DO i = 1, 10000
  READ(11,*,END=3) spex(i),sper(i)
  nspe = nspe + 1
END DO
3 CONTINUE
CLOSE(11)

!! WRITE(*,*) 'SPE file lines:', nspe

! Test if speresp function works properly
!DO i = 1, 50000
!  spet = speresp(1.0,spex,sper,nspe)
!  WRITE(12,*) spet
!END DO
!CLOSE(12)

! Treat input Cherenkov photon number as a rate per simulated interval
! and randomize the number of arriving photons according to Poisson distribution

IF (ipor == 1) THEN

  racer = REAL(ncer)
  DO i = 1, 100000
    CALL pocher(racer,deltime,ncerpo)
    pocer(i) = ncerpo
  END DO

  pomin = MINVAL(pocer(1:100000))
  pomax = MAXVAL(pocer(1:100000))

  ponu = 0
  DO i = 1, 100000
    DO j = pomin,pomax
      IF (pocer(i) == j) ponu(j) = ponu(j)+1
    END DO
  END DO

  amsum = 0.0
  DO j = pomin, pomax
! random.photons
!    WRITE(55,*) j, ponu(j)
    amsum = amsum + REAL(ponu(j)*j)
  END DO 
!  CLOSE(55)
  ramean = amsum/100000.0

  amsum = 0.0
  DO j = pomin, pomax
    amsum = amsum + (ramean-REAL(j))*(ramean-REAL(j))*REAL(ponu(j))/100000.0
  END DO
  ravar = amsum
  rasigma = SQRT(amsum)

  WRITE(*,*) 'Randomized Cherenkov photon number', ncerpo, ramean, rasigma, ravar

  ncer = ncerpo

END IF

IF (itim == 1) THEN
! Generate random arrival times of photons
! Add jitter time
  DO i = 1, ncer
    atime(i) = bodi(time0,deltime)
!  WRITE(*,*) atime(i)
    atime(i) = atime(i) + rano(0.0,tjitter)
!!  WRITE(*,*) atime(i)
  END DO 
END IF

IF (itim == 0) THEN
! Instrumental resolution tests, all photons arrive at t = 100 ns
  DO i = 1, ncer
    atime(i) = time0
  END DO
END IF

atimemin = MINVAL(atime(1:ncer))
atimemax = MAXVAL(atime(1:ncer))
IF (ipr == 1) WRITE(*,*) 'Time range:', atimemin,atimemax

sumati = 0.0
DO i = 1, ncer
  atime(i) = atime(i) - atimemin
  sumati = sumati + atime(i)
!  WRITE(*,*) i, atime(i)
END DO

artime = sumati/REAL(ncer)
detime = atimemax - atimemin

sumati = 0.0
DO i = 1, ncer
  sumati = sumati + (atime(i)-artime)**2.0
END DO

sitime = SQRT(sumati/REAL(ncer))

WRITE(10,*) 'Photons:', ncer, artime, sitime, detime

DO i = 1, ncer
  cerart(i) = atime(i) + atimemin
END DO

IF (ipr == 1) THEN
  DO i = 1, ncer
! cherenkov.photons
    WRITE(21,*) cerart(i)
  END DO
  CLOSE(21)
END IF

! Simulate NSB photons arrival times

CALL nsbat(nsbrate,tstart,tstop,nsbnum,nsbart)

CALL tisort(nsbart,nsbnum)

IF (ipr == 1) WRITE(*,*) 'NSB photons:', nsbnum

!nsbnum = 0
!DO i = 1, 10000
!  nsbnum = nsbnum + 1
!  READ(22,*,END=10) nsbart(i)
!END DO
!10 CONTINUE

IF (ipr == 1) THEN
  DO i = 1, nsbnum
! nsb.photons
    WRITE(22,*) nsbart(i)
!!  WRITE(*,*) i, nsbart(i)
  END DO
  CLOSE(22)
END IF

! Simulate thermal noise pulses arrival times

CALL nsbat(tenrate,tstart,tstop,tennum,tenart)

CALL tisort(tenart,tennum)

IF (ipr == 1) WRITE(*,*) 'Noise pulses:', tennum

!nsbnum = 0
!DO i = 1, 10000
!  nsbnum = nsbnum + 1
!  READ(22,*,END=10) nsbart(i)
!END DO
!10 CONTINUE

IF (ipr == 1) THEN
  DO i = 1, tennum
! thermal.photons
    WRITE(23,*) tenart(i)
  END DO
  CLOSE(23)
END IF

! Merge Cherenkov and NSB photons into photon list
! DIFFERENT QUANTUM EFFICIENCY (DUE TO DIFFERENT SPECTRA) CAN BE TAKEN INTO ACCOUNT
np = 1
DO i = 1, nsbnum
  ptime(np) = nsbart(i)
  np = np + 1
END DO

DO i = 1, ncer
  ptime(np) = cerart(i)
  np = np + 1
END DO

DO i = 1, tennum
  ptime(np) = tenart(i)
  np = np + 1
END DO
np = np - 1

CALL tisort(ptime,np)

IF (ipr == 1) WRITE(*,*) 'Total photons:', np

!! DO i = 1, np
!!   WRITE(*,*) i, ptime(i)
!! END DO

! Investigated time period
! Finer binning to avoid 0.1 ns grouping of photon pulses
! 0.01 ns instead, well below time resolution of GAPD

!nt = 2000
DO i = 1, nt
  tt(i) = tstart + (REAL(i-1)+0.5)*dt
  pude(i) = 0.0
END DO

DO i = 1, nfine
  t1(i) = tstart + REAL(i-1)*dtfine
  t2(i) = t1(i) + dtfine
  tfine(i) = t1(i) + 0.5*dtfine
  nf(i) = 0
  pudefine(i) = 0.0
END DO

! Sum photons over time bins

!DO i = 1, nt
DO i = 1, nfine
  DO j = 1, np
    IF (ptime(j) >= t1(i) .AND. ptime(j) < t2(i)) nf(i) = nf(i) + 1
  END DO
END DO

ne = 0
tfines = tstop
ifines = 1
DO i = 1, nfine
  IF (nf(i) > 0) THEN
! detector.signals
IF (ipr == 1) WRITE(24,*) i, tfine(i), nf(i)
    ne = ne + 1
    te(ne) = tfine(i)
    IF (ne == 1) THEN
      ifines = i
      tfines = te(ne)
    END IF
    me(ne) = nf(i)
  END IF
END DO
IF (ipr == 1) CLOSE(24)

!WRITE(*,*) 'tfines,ifines:', tfines,ifines

! Signal amplitude blurred according to single photoelectron spectrum
! Each photon amplitude blurred separately and then summed up
DO j = 1, ne
  am(j) = 0.0
  DO k = 1, me(j)
IF (ispe == 1) am(j) = am(j) + speresp(onepe,spex,sper,nspe)
IF (ispe == 0) am(j) = REAL(me(j))*onepe
  END DO
! photon.bunches
IF (ipr == 1) WRITE(53,*) te(j),me(j),am(j)
END DO
IF (ipr == 1) CLOSE(53)

! Sum elementary pulses into detector pulse ('signal' in figure on page 37)

! 250 MHz
sh = 1.4154
uh = 0.83725
dth = 2.95
! Silicon detectors pulse shape 
! S10985-050C
!tra = 1.12
!tde = 25.33
! S12516-050C
tra = 3.38
tde = 52.53
! test
!tde = 30.0
! S12516-100C
!tra = 2.35
!tde = 92.80
dth = 0.0
! 1 GHz
!sh = 1.425
!uh = 0.90954
!dth = 3.01
DO i = ifines-1, nfine
  DO j = 1, ne
    IF (te(j) <= tfine(i)) THEN
      th = te(j) + dth
      ah = REAL(me(j))*onepe
! PMT signal
!      pude(i) = pude(i) + hyga(tt(i),ah,th,sh,uh)
! GAPD signal
      pudefine(i) = pudefine(i) + dexe(tfine(i),am(j),th,tra,tde)
    END IF
  END DO
END DO

DO i = ifines-1,ifines+1000
  WRITE(52,*) tfine(i),pudefine(i)
END DO
CLOSE(52)

CALL curin(nfine,nt,tfine,tt,pudefine,pude)

desum = 0.0
ide = 0
DO i = 1, nt
  desum = desum + pude(i)
  IF (desum < 0.0001 ) ide = i
! detector.pulse
IF (ipr == 1) WRITE(31,*) tt(i), pude(i)
END DO
IF (ipr == 1) CLOSE(31)

! Reference pulse integral
resum = desum*dt
pudemax = MAXVAL(pude)

WRITE(10,*) 'Detector pulse integral and maximum:', resum, pudemax
IF (ipr == 1) WRITE(*,*) 'Detector pulse integral and maximum:', resum, pudemax, ide

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
pufi = 0.0
pufine = 0.0
! Only for tamp computation
!DO i = ifines-20, ifines+30
!  pufine(i) = amfire(i,pudefine,tfine,af1,af2,af3)
!  IF (pufine(i) > 0.001 .AND. pufine(i-1) < 0.001) tamp = tfine(i)
!  WRITE(51,*) tfine(i),pufine(i)
!END DO
!CLOSE(51)

! S12516, both 50C and 100C
tamp = 900.015
!tamp = 900.025

DO i = ide, nt
  pufi(i) = amfire(i,pude,tt,af1,af2,af3)
  fisum = fisum + pufi(i)
END DO

IF (ipr == 1) WRITE(*,*) 'Amplifier signal starts at:', tamp

! Photons arrive to ADC at random ADC phase
! Sampling starts at tstart + 0.75*tsam (see sampling algorithm)
! Therefore ADC start time is in the range (0.75*tsam,1.75*tsam)
! Random phase
IF (ipha == 1) tadc = bodi(tstart+1.25*tsam,tsam)
! Fixed phase for amplitude
!IF (ipha == 0) tadc = tstart + 1.243*tsam
! Fixed phase for time
IF (ipha == 0) tadc = tstart + 1.502*tsam
na = 1
tad = tadc
DO WHILE (tad < tt(nt)-0.5*tsam)
  ta(na) = tadc + tsam*(REAL(na)-0.5)
  tad = ta(na)
  na = na + 1
END DO
na = na - 1
IF (ipr == 1) WRITE(*,*) 'ADC start time, na:', tadc, na, tad

! Fine interpolation bins
ni = 1
tic = tadc
DO WHILE (tic < tad)
  ti(ni) = tadc + REAL(ni)*dti
  tic = ti(ni)
  ni = ni + 1
END DO
ni = ni-1
IF (ipr == 1) WRITE(*,*) 'Fine interpolation bins ni:', ni,ti(1),tic

! Finest interpolation bins
nfine = 1
tic = tadc
ifine0 = 1
DO WHILE (tic < tad)
  tfine(nfine) = tadc + REAL(nfine)*dtfine
  tic = tfine(nfine)
  IF (tic < time0) ifine0 = nfine
  nfine = nfine + 1
END DO
nfine = nfine-1
! Index range for MAXVAL search, used for deconvolved and shaped signals
! time0-10 ns -- time0+30 ns
ifine1 = ifine0 - INT(10.0/dtfine)
ifine2 = ifine0 + INT(30.0/dtfine)
!WRITE(*,*) 'Finest interpolation bins nfine:', nfine,tfine(1),tic
!WRITE(*,*) 'Finest interpolation central ranges:', ifine1,tfine(ifine1),ifine0,tfine(ifine0),ifine2,tfine(ifine2)

! white.noise added to amplified and filtered signal
amsum = 0.0
fasum = 0.0
punor = desum/fisum
DO i = 1, nt
  wnoisa = rano(wnom,wnos)
! white.noise
!  WRITE(54,*) tt(i),wnoisa
  pufi(i) = punor*pufi(i)
  pufi(i) = pufi(i) + wnoisa
  amsum = amsum + wnoisa*wnoisa
  fasum = fasum + 1
END DO
!CLOSE(54)
wnomean = SQRT(amsum/fasum)
wnomeanpe = wnomean/onepa
IF (ipr == 1) WRITE(*,*) 'White noise amplitude, absolute and in p.e. units:', wnomean, wnomeanpe

IF (ipr == 1) THEN
  DO i = 1, nt
!  pufi(i) = desum*pufi(i)/fisum
! amplifier.pulse
    WRITE(32,*) tt(i), pufi(i)
  END DO
  CLOSE(32)
END IF
pufisum = SUM(pufi)*dt
pufimax = MAXVAL(pufi(1:nt))

WRITE(10,*) 'Amplifier pulse integral and maximum:', pufisum, pufimax

! Amplitude and time at signal maximum
imax = MAXLOC(pufi(1:nt))
maxtim = tt(imax(1))

WRITE(10,*) 'Amplifier signal, arrival time:', maxtim

IF (ipr == 1) WRITE(*,*) 'Amplifier pulse integral and maximum:', pufisum, pufimax

! Finest interpolation before digitization
CALL curin(nt,nfine,tt,tfine,pufi,afine)

!DO i = 30000,50000
!  WRITE(51,*) tfine(i),afine(i),i
!END DO
!CLOSE(51)

! Reference time for ADC phase, when the ADC input signal (amplifier pulse) starts

! FADC sampling (SAMPLING, NOT REBINNING, SOME INFORMATION LOST)
! Average over a half clock cycle, starting 0.75 cycle before ta(i)
! because average is done over past pulse (-0.25,+0.25)
! at the ADC time bin start
! Rounding to the nearest integer
fadc = 0.0
DO i = 1, na
  j1 = ANINT((ta(i)-tsam34-tfine(1))/dtfine)+1
  IF (j1 < 1) j1 = 1
  IF (j1 > nfine) j1 = nfine
  j2 = ANINT((ta(i)-tsam14-tfine(1))/dtfine)
  IF (j2 > nfine) j2 = nfine
  tam1 = ta(i)-tsam34
  tam2 = ta(i)-tsam14
  amsum = 0.0
  DO j=j1,j2
    amsum = amsum + afine(j)
  END DO
  fasum = amsum/REAL(j2-j1+1)
!  fadc(i) = REAL(ANINT(fasum))
!  IF (ipr == 1) WRITE(44,*) ta(i),fadc(i)  
  WRITE(*,*) 'Fasum 1', fasum
! FADC noise according to sim_telarray
  IF (fanoise > 0.0) fasum = rano(fasum,fanoise)
  WRITE(*,*) 'Fasum 2', fasum
! Pedestal variations affecting signal (e.g. due to temperature variations)
  IF (pevar > 0.0) fasum = rano(fasum,pevar)
  WRITE(*,*) 'Fasum 3', fasum
  fadc(i) = REAL(ANINT(fasum))
! fadc.samples
IF (ipr == 1) WRITE(33,*) ta(i),fadc(i)  
END DO
IF (ipr == 1) CLOSE(33)

! FADC samples time is SHIFTED AHEAD by 0.5 tsam
! Transform to the real time for a further signal analysis
sampsum = 0.0
DO i = 1, na
  ta(i) = ta(i) - tsam12
  sampsum = sampsum + fadc(i)*tsam
END DO
sampmax = MAXVAL(fadc(1:na))

WRITE(10,*) 'FADC samples integral and maximum:', sampsum, sampmax
IF (ipr == 1) WRITE(*,*) 'FADC samples integral and maximum:', sampsum, sampmax

! Smoothing the ADC counts with a moving average of 1 clock cycle

! First linear interpolation is done to obtain ADC amplitude for 1 ns bins
! Linear and parabolic interpolation of the sampled signal

fadi = 0.0
CALL curin(na,ni,ta,ti,fadc,fadi)

isamsum = SUM(fadi)*dti
isammax = MAXVAL(fadi(1:ni))

WRITE(10,*) 'Linear interpolation of FADC samples integral and maximum:', isamsum, isammax

! Moving average of interpolated ADC signal

fada = 0.0
CALL movave(ni,nt,ti,tt,fadi,fada,tsam)

misasum = 0.0
DO i = 1, nt
  misasum = misasum + ABS(fada(i))
! fadc.aver
IF (ipr == 1) WRITE(35,*) tt(i),fada(i)
END DO
IF (ipr == 1) CLOSE(35)

misasum = misasum*dt
misamax = MAXVAL(fada(1:nt))

WRITE(10,*) 'Moving average of interpolated FADC samples integral and maximum:', misasum, misamax
 
! Smoothing the averaged FADC signal with another interpolation

DO i = 1, nt
  tu(i) = tt(i) - 0.5*dt
END DO

fadu = 0.0
CALL curin(nt,nt,tt,tu,fada,fadu)

! Derivative of the smoothed averaged FADC signal
fade = 0.0
dt12 = 1.0/(12.0*dt)
DO i = 1, nt
   IF (i > 2 .AND. i < nt-1) THEN
     fade(i) = (8.0*fada(i+1)+fada(i-2)-8.0*fada(i-1)-fada(i+2))*dt12
   END IF
END DO
     
IF (ipr == 1) THEN
  DO i = 1, nt
! fadc.dave
    WRITE(37,*) tt(i),fade(i)
  END DO
  CLOSE(37)      
END IF

! Fine interpolation of the averaged signal 

fadar = 0.0
CALL curin(nt,ni,tt,ti,fade,fadi)
CALL movave(ni,nt,ti,tt,fadi,fadar,tsam)

amsum = 0.0
DO i = 1, nt
! fadc.daver
  amsum = amsum + ABS(fadar(i))
  fadar(i) = defa*signorm*fadar(i)
IF (ipr == 1) WRITE(39,*) tt(i),fadar(i)
END DO
IF (ipr == 1) CLOSE(39)
amsum = amsum*dt

IF (ipr == 1) WRITE(*,*) 'Derivative integral:', amsum
WRITE(10,*) 'Derivative integral:', amsum

! Deconvolved signal: sum of averaged signal and averaged derivative
fadeco = 0.0
DO i = 2, nt
  fadeco(i) = (fadu(i-1)+fadu(i))/2.0 + fadar(i)
! fadc.deco
IF (ipr == 1) WRITE(41,*) tt(i),fadeco(i)
END DO
IF (ipr == 1) CLOSE(41)

! Amplitude and time at signal maximum
!maxamp = MAXVAL(fadeco(1:nt))
!imax = MAXLOC(fadeco(1:nt))
!maxtim = tt(imax(1))

!IF (ipr == 1) WRITE(*,*) 'Deconvolved signal maximum and arrival time, coarse:', maxamp,maxtim

! Fine interpolation for a better estimate of maximum amplitude and time
CALL curin(nt,nfine,tt,tfine,fadeco,afine)
imax = MAXLOC(afine(ifine1:ifine2))   
maxamp = afine(imax(1)+ifine1-1)
maxtim = tfine(imax(1)+ifine1-1)

IF (ipr == 1) WRITE(*,*) 'Deconvolved signal maximum and arrival time, fine:', maxamp,maxtim
  
sisum = SUM(fadeco)*dt

! Background level

IF (pedelevel < 0.0) THEN
  bgsumde = 0.0
  nbgde = 0
  DO i = 1, nt
    IF (tt(i) > bi1 .AND. tt (i) < bo1) THEN
      bgsumde = bgsumde + fadeco(i)
      nbgde = nbgde + 1
    END IF
    IF (tt(i) > bi2 .AND. tt (i) < bo2) THEN
      bgsumde = bgsumde + fadeco(i)
      nbgde = nbgde + 1
    END IF
  END DO
  bgsumde = bgsumde/REAL(nbgde)
ELSE
  bgsumde = pedelevel  
END IF

! Additional uncertainty parameters similar to those included in sim_telarray
! Simplified approach (we do not take into account that some parameters in 
! sim_telarray appear twice, in both additive (pedestal) and multipliticative 
! (calibration) factors

! Pedestal parameters
! Channel to channel variations
WRITE(*,*) 'Pedestal 0', bgsumde
IF (pepar1 > 0.0) bgsumde = rano(bgsumde,pepar1)
WRITE(*,*) 'Pedestal 1', bgsumde
! FADC deviations for a single channel
IF (pepar2 > 0.0) bgsumde = rano(bgsumde,pepar2)
WRITE(*,*) 'Pedestal 2', bgsumde
! Error of initial calibration
IF (pepar3 > 0.0) bgsumde = rano(bgsumde,pepar3)
WRITE(*,*) 'Pedestal 3', bgsumde

! Calibration (amplitude to npe) parameters
! Quantum efficiency variations
WRITE(*,*) 'Calibration 0', refeped
IF (qevar > 0.0) refeped = rano(refeped,qevar*refeped)
WRITE(*,*) 'Calibration 1', refeped
! Gain variations
IF (gavar > 0.0) refeped = rano(refeped,gavar*refeped)
WRITE(*,*) 'Calibration 2', refeped
! Sensitivity variations
IF (sevar > 0.0) refeped = rano(refeped,sevar*refeped)
WRITE(*,*) 'Calibration 3', refeped

WRITE(10,*) 'Background level (deconvolved):', bgsumde
IF (ipr == 1) WRITE(*,*) 'Background level (deconvolved):', bgsumde

! Reconstructed number of photoelectrons
repe = (maxamp-bgsumde)/refeped

IF (ipr == 1) WRITE(*,*) 'Deconvolved signal, amplitude:', maxamp, sisum, sampsum, resum
IF (ipr == 1) WRITE(*,*) 'Deconvolved signal, arrival time:', maxtim,ptime(1)
IF (ipr == 1) WRITE(*,*) 'Deconvolved signal, number of photoelectrons:', repe

WRITE(10,*) 'Deconvolved signal, amplitude:', maxamp, sisum, sampsum, resum
WRITE(10,*) 'Deconvolved signal, arrival time:', maxtim,ptime(1)
WRITE(10,*) 'Deconvolved signal, number of photoelectrons:', repe

! Shaped signal: three moving averages of the deconvolved signal

! 250 MHz
tave = 2.0*tsam
! 1 GHz
!tave = 8.0*tsam
!tave = 8.0

fasa = 0.0
CALL curin(nt,ni,tt,ti,fadeco,fadin)
CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

CALL curin(nt,nt,tt,tu,fasa,fadu)
CALL movave(nt,nt,tu,tt,fadu,fasa,tave)

CALL curin(nt,nt,tt,tu,fasa,fadu)
CALL movave(nt,nt,tu,tt,fadu,fasa,tave)

CALL curin(nt,nt,tt,tu,fasa,fadu)
CALL movave(nt,nt,tu,tt,fadu,fasa,tave)

! Two more for 1 GHz
!CALL herin(nt,ni,tt,ti,fasa,fadin)
!CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

!CALL herin(nt,ni,tt,ti,fasa,fadin)
!CALL movave(ni,nt,ti,tt,fadin,fasa,tave)

!maxamp = MAXVAL(fasa(1:nt))
!imax = MAXLOC(fasa(1:nt))
!maxtim = tt(imax(1))

!IF (ipr == 1) WRITE(*,*) 'Shaped signal maximum and arrival time, coarse:', maxamp,maxtim

! Fine interpolation for a better estimate of maximum amplitude and time
CALL curin(nt,nfine,tt,tfine,fasa,afine)
imax = MAXLOC(afine(ifine1:ifine2))
maxtim = tfine(imax(1)+ifine1-1)   
maxamp = afine(imax(1)+ifine1-1)   

IF (ipr == 1) WRITE(*,*) 'Shaped signal maximum and arrival time, fine:', maxamp,maxtim

sisum = SUM(fasa)*dt

IF (ipr == 1) THEN
  DO i = 1, nt
! fadc.shaped
    WRITE(42,*) tt(i),fasa(i)
  END DO
  CLOSE(42)
END IF

! Background level
IF (pedelevel < 0.0) THEN
  bgsumsh = 0.0
  nbgsh = 0
  DO i = 1, nt
    IF (tt(i) > bi1 .AND. tt (i) < bo1) THEN
      bgsumsh = bgsumsh + fasa(i)
      nbgsh = nbgsh + 1
    END IF
    IF (tt(i) > bi2 .AND. tt (i) < bo2) THEN
      bgsumsh = bgsumsh + fasa(i)
      nbgsh = nbgsh + 1
    END IF
  END DO
  bgsumsh = bgsumsh/REAL(nbgsh)
ELSE
  bgsumsh = pedelevel
END IF

! Additional uncertainty parameters similar to those included in sim_telarray
! Simplified approach (we do not take into account that some parameters in 
! sim_telarray appear twice, in both additive (pedestal) and multipliticative 
! (calibration) factors

! Pedestal parameters
! Channel to channel variations
!WRITE(*,*) 'Pedestal 0', bgsumde
IF (pepar1 > 0.0) bgsumsh = rano(bgsumsh,pepar1)
!WRITE(*,*) 'Pedestal 1', bgsumde
! FADC deviations for a single channel
IF (pepar2 > 0.0) bgsumsh = rano(bgsumsh,pepar2)
!WRITE(*,*) 'Pedestal 2', bgsumde
! Error of initial calibration
IF (pepar3 > 0.0) bgsumsh = rano(bgsumsh,pepar3)
!WRITE(*,*) 'Pedestal 3', bgsumde

! Calibration (amplitude to npe) parameters
! Quantum efficiency variations
!WRITE(*,*) 'Calibration 0', refeped
IF (qevar > 0.0) refepes = rano(refepes,qevar*refepes)
!WRITE(*,*) 'Calibration 1', refeped
! Gain variations
IF (gavar > 0.0) refepes = rano(refepes,gavar*refepes)
!WRITE(*,*) 'Calibration 2', refeped
! Sensitivity variations
IF (sevar > 0.0) refepes = rano(refepes,sevar*refepes)
!WRITE(*,*) 'Calibration 3', refeped

WRITE(10,*) 'Background level (shaped):', bgsumsh
IF (ipr == 1) WRITE(*,*) 'Background level (shaped):', bgsumsh

! Reconstructed number of photoelectrons
repe = (maxamp-bgsumsh)/refepes

IF (ipr == 1) WRITE(*,*) 'Shaped signal, amplitude:', maxamp, sisum, sampsum, resum
IF (ipr == 1) WRITE(*,*) 'Shaped signal, arrival time:', maxtim,ptime(1)
IF (ipr == 1) WRITE(*,*) 'Shaped signal, number of photoelectrons:', repe

WRITE(10,*) 'Shaped signal, amplitude:', maxamp, sisum, sampsum, resum
WRITE(10,*) 'Shaped signal, arrival time:', maxtim,ptime(1)
WRITE(10,*) 'Shaped signal, number of photoelectrons:', repe

! Phase difference
! 0 when the end of sampled interval (-0.25,+0.25) is at time 
! when the ADC input signal becomes non-zero
! 1 (tsam) when the previous sampled interval approaches non-zero input
! phase = ta(i)-0.25tsam - tinput0
! tinput0 different from Cherenkov photon arrival time time0 because
! the true sampled signal is that after amplifier and filter, thus
! delayed, also detector should give some delay, even silicon detector
! for which it is assumed to be 0
!nbgsh = FLOOR((time0-tadc+0.75*tsam)/tsam)+1
nbgsh = FLOOR((tamp-tadc+0.75*tsam)/tsam)+1
sisum = ta(nbgsh)+0.5*tsam
maxamp = sisum-0.75*tsam
maxtim = sisum-0.25*tsam
phadc = maxtim-tamp
WRITE(10,*) 'FADC delay time and phase:', tadc,phadc,maxamp,maxtim,tamp
IF (ipr == 1) WRITE(*,*) 'FADC delay time and phase:', tadc,nbgsh,ta(nbgsh),sisum,maxamp,maxtim,phadc

call ETIME(runtim2,telaps2)

ietime = runtim2(1)
      
IF (ipr == 1) WRITE(*,400) ietime

400 FORMAT(1x,'Execution time:',i9,' seconds')

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

! Cherenkov photons random Poisson number
SUBROUTINE pocher(nce,tidel,npoch)

REAL, INTENT(IN) :: nce,tidel
INTEGER, INTENT(OUT) :: npoch

REAL :: raa,rab,lambda2,rate
REAL(KIND=8) :: lambda,nsbnur,wran,pran
INTEGER :: i

! Find Poisson distributed number of Cherenkov photons
! Cherenkov photon rate converted to photons per ns
rate = nce
lambda = REAL(rate,KIND=8)
lambda2 = 6.0*rate
DO
  CALL RANDOM_NUMBER(raa)
  npoch = INT(lambda2*raa)
  nsbnur = REAL(npoch,KIND=8)
  wran = dpois(nsbnur,lambda)
  pran = DEXP(wran)
  CALL RANDOM_NUMBER(rab)
  IF (rab < pran) EXIT
END DO

RETURN

END SUBROUTINE pocher

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
REAL(KIND=8) :: denx

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
    denx = 1.0/(xi(j3)-xi(j2))
    xo2 = (xo(i)-xi(j2))*(xo(i)-xi(j2))*denx*denx
    xo3 = (xo(i)-xi(j2))*(xo(i)-xi(j2))*(xo(i)-xi(j2))*denx*denx*denx
!    yy1 = (fi(j2)-fi(j1))*(1.0+bias)*(1.0-tension)/2.0
!    yy1 = yy1 + (fi(j3)-fi(j2))*(1.0-bias)*(1.0-tension)/2.0
!    yy2 = (fi(j3)-fi(j2))*(1.0+bias)*(1.0-tension)/2.0
!    yy2 = yy2 + (fi(j4)-fi(j3))*(1.0-bias)*(1.0-tension)/2.0
    yy1 = 0.5*(fi(j2)-fi(j1))
    yy1 = yy1 + 0.5*(fi(j3)-fi(j2))
    yy2 = 0.5*(fi(j3)-fi(j2))
    yy2 = yy2 + 0.5*(fi(j4)-fi(j3))
    x1 = 2.0*xo3 - 3.0*xo2 + 1.0
    x2 = xo3 - 2.0*xo2 + (xo(i)-xi(j2))*denx
    x3 = xo3 - xo2
    x4 = -2.0*xo3 + 3.0*xo2
    fo(i) = x1*fi(j2) + x2*yy1 + x3*yy2 + x4*fi(j3)
!    IF (xo(i) > 105.0 .AND. xo(i) < 106.0) THEN
!      WRITE(*,*) 'herin',xo(i),yy1,yy2,x1,x2,x3,x4,fo(i)
!    END IF
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
REAL :: x1,x2,x3,x4,y1,y2,y3,y4,dxi,dxo,xo2,xo3,yy1,yy2,denx

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,j4

DO i = 1, no
  fo(i) = 0.0
  IF ( xo(i) .GE. xi(2) .AND. xo(i) .LE. xi(ni-1) ) THEN
    DO j = 2, ni-1
      IF (xi(j) <= xo(i) .AND. xi(j+1) > xo(i)) THEN
        j1 = j-1
        j2 = j
        j3 = j+1
        j4 = j+2
      END IF
    END DO
    denx = (xo(i)-xi(j2))/(xi(j3)-xi(j2))
    xo2 = denx*denx
    xo3 = xo2*denx
    x1 = -0.5*fi(j1) + 1.5*fi(j2) - 1.5*fi(j3) + 0.5*fi(j4)
    x2 = fi(j1) - 2.5*fi(j2) + 2.0*fi(j3) - 0.5*fi(j4)
    x3 = -0.5*fi(j1) + 0.5*fi(j3)
    x4 = fi(j2)
    fo(i) = x1*xo3 + x2*xo2 + x3*denx + x4
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
REAL :: asum,wsum,dfo,dxi,we1,we2,we3,we4,dts,dexi

INTEGER, INTENT(IN) :: ni,no
INTEGER :: i,j,j1,j2,j3,dj

dts = 0.5*ts
dxi = 0.5*(xi(ni)-xi(1))/REAL(ni-1)
dexi = 0.5/dxi

DO i = 1, no
  fo(i) = 0.0
  IF (xo(i)-dts > xi(1) .AND. xo(i)+dts < xi(ni)) THEN
    DO j = 1, ni
      IF (xo(i)-dts > xi(j) .AND. xo(i)-dts <= xi(j+1)) j1 = j+1
      IF (xo(i)+dts >= xi(j) .AND. xo(i)+dts < xi(j+1)) j2 = j
    END DO
    asum = 0.0
    wsum = 0.0
    we1 = 0.0
    we2 = 0.0
    we3 = 0.0
    we4 = 0.0
! Partial overlap 
    IF (xo(i)-dts < xi(j1-1)+dxi) THEN
      we1 = (xi(j1-1)+dxi-xo(i)+dts)*dexi
      asum = asum + we1*fi(j1-1) + fi(j1)
      wsum = wsum + we1 + 1.0
    END IF
    IF (xo(i)-dts >= xi(j1-1)+dxi) THEN
      we2 = (xi(j1)+dxi-xo(i)+dts)*dexi
      asum = asum + we2*fi(j1)
      wsum = wsum + we2
    END IF
    DO j = j1+1, j2-1
      asum = asum + fi(j)
      wsum = wsum + 1.0
    END DO
    IF (xo(i)+dts > xi(j2+1)-dxi) THEN
      we3 = (xi(j2+1)+dxi-xo(i)-dts)*dexi
      asum = asum + we3*fi(j2+1) + fi(j2)
      wsum = wsum + we3 + 1.0
    END IF
    IF (xo(i)+dts <= xi(j2)+dxi) THEN
      we4 = (xi(j2)+dxi-xo(i)-dts)*dexi
      asum = asum + we4*fi(j2)
      wsum = wsum + we4
    END IF
!    dj = j2-j1+1
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
  REAL :: speresp,ra,ram,rasp,dspx
  INTEGER :: ira
  dspx = 1.0/(spx(2)-spx(1))
  DO
    CALL RANDOM_NUMBER(ra)
    ram = spx(1) + ra*(spx(ns)-spx(1))
    ira = ANINT((ram-spx(1))*dspx)
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
  REAL :: rano,nra,dino,sig2,nos12,rando
  sig2 = 1.0/(2.0*nos*nos)
  nos12 = 12.0*nos
  DO
    CALL RANDOM_NUMBER(nra)
    rando = nos12*(nra-0.5)
    dino = EXP(-rando*rando*sig2)
    CALL RANDOM_NUMBER(nra)
    IF (nra < dino) EXIT
  END DO
  rano = rando + noc
  RETURN
END FUNCTION rano

! Convolution with the amplifier/filters response
FUNCTION amfire(ipu,pud,tid,ama,amb,amt)
  REAL, INTENT(IN) :: ama,amb,amt
  REAL, DIMENSION(200000), INTENT(IN) :: pud,tid
  REAL :: amfire,ares,sumco,detu,dett,ramt
  INTEGER, INTENT(IN) :: ipu
  INTEGER :: i
  sumco = 0.0
  dett = tid(2)-tid(1)
  ramt = 1.0/amt
  DO i = 1, ipu
    detu = tid(i)-tid(1)
    ares = (ama+amb*detu)*EXP(-detu*ramt)
    sumco = sumco + pud(ipu-i)*ares
  END DO
  amfire = sumco*dett
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


END PROGRAM cadisi
