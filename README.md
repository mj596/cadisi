cadisi
======
== Cadisi Notes ==

----

* gfortran cadisi.f90 -o cadisi -g
* Cadisi needs 'cadisi.par' in main catalogue to run. Else: segfault.
* Cadisi needs 'mppcspenormean.dat' in main catalogue to run. Else: segfault.
* 'mppcspenormean.dat' - normalized photodetector SPE (first max not at 1 pe!)
* Results file: 'cadisi.res'
* No other files are opened nor requested by cadisi:
  {{{
  open("cadisi.par", O_RDWR|O_CREAT, 0666) = 3
  open("cadisi.results", O_RDWR|O_CREAT, 0666) = 3
  open("mppcspenormean.dat", O_RDWR|O_CREAT, 0666) = 4
  }}}
* For
   '1  ipr: Printing flag'
  cadisi opens and writes:
  {{{
  open("cadisi.par", O_RDWR|O_CREAT, 0666) = 3
  open("cadisi.results", O_RDWR|O_CREAT, 0666) = 3
  open("mppcspenormean.dat", O_RDWR|O_CREAT, 0666) = 4
  open("cherenkov.photons", O_RDWR|O_CREAT, 0666) = 5
  open("nsb.photons", O_RDWR|O_CREAT, 0666) = 6
  open("thermal.photons", O_RDWR|O_CREAT, 0666) = 7
  open("detector.signals", O_RDWR|O_CREAT, 0666) = 8
  open("detector.pulse", O_RDWR|O_CREAT, 0666) = 9
  open("amplifier.pulse", O_RDWR|O_CREAT, 0666) = 10
  open("fadc.samples", O_RDWR|O_CREAT, 0666) = 11
  open("fadc.aver", O_RDWR|O_CREAT, 0666) = 12
  open("fadc.dave", O_RDWR|O_CREAT, 0666) = 13
  open("fadc.daver", O_RDWR|O_CREAT, 0666) = 14
  open("fadc.deco", O_RDWR|O_CREAT, 0666) = 15
  open("fadc.shaped", O_RDWR|O_CREAT, 0666) = 16
  open("fadc.shapin", O_RDWR|O_CREAT, 0666) = 17
  open("flasim.test", O_RDWR|O_CREAT, 0666) = 18
  open("photon.bunches", O_RDWR|O_CREAT, 0666) = 19
  }}}
  * Created files:
  cherenkov.photons - ?
  nsb.photons - ?
  thermal.photons - ?
  Possibly number of photons vs arrival time ?

  detector.signals - ?
  detector.pulse - possibly detector amplitude (mV?) vs time (ns?). Makes sense as in cadisi.par:
  {{{
	0.0          tstart: Start time of computed period (ns)
	2000.0       tstop: Stop time of computed period (ns)
  }}}
  and time axis spans 0-2000 range.

  amplifier.pulse - possibly amplifier amplitude (mV?) vs time (ns?). It looks like photodetector pulse shifted by a couple of ns back in time and with noise added.

  fadc.samples - possibly digitized signal from pre-amplifier
  fadc.aver - averaged (somehow) digitied signal
  fadc.dave - a derivative of averaged digitized signal
  fadc.daver - looks like a digitized signal derivative but vertically streched to match (something??)
  fadc.deco - 'deco' suggests deconvolution
  fadc.shaped - well... shaped signal
  fadc.shapin - NO DATA
  flasim.test - NO DATA
  photon.bunches - ??



----
