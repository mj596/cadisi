#!/usr/bin/perl

# FlasCam simulation tools

use POSIX;

$ni = $ARGV[0];

open (Difile,'>','flaszin.dis');

#$ni = 3000;

$suma = 0;
$sumt = 0;
$minram = 1000;
$maxram = -1000;
$minti = 1000;
$maxti = -1000;
for ($i=0;$i<$ni;$i++) {
  system "./cadisint";
  open (Refile,'<','cadisint.results');
  while (<Refile>) {
    chomp;
    @Line = split;
    if ($Line[0] eq "Deconvolved" && $Line[2] eq "amplitude:") {$dam[$i] = $Line[3]}
    if ($Line[0] eq "Shaped" && $Line[2] eq "amplitude:") {$ram[$i] = $Line[3]}
    if ($Line[0] eq "Deconvolved" && $Line[2] eq "arrival") {$tide[$i] = $Line[4]}
    if ($Line[0] eq "Shaped" && $Line[2] eq "arrival") {$tish[$i] = $Line[4]}
    if ($Line[0] eq "Deconvolved" && $Line[2] eq "fine") {
      $rame[$i] = $Line[7];
      $tise[$i] = $Line[8];
    }
    if ($Line[0] eq "Shaped" && $Line[2] eq "fine") {
      $rami[$i] = $Line[7];
      $tisi[$i] = $Line[8];
    }
    if ($Line[0] eq "Photons:") {
      $tar = $Line[2];
      $sar = $Line[3];
      $dar = $Line[4];
    }
    if ($Line[0] eq "Background" && $Line[2] eq "(deconvolved):") {$badec[$i] = $Line[3]}
    if ($Line[0] eq "Background" && $Line[2] eq "(shaped):") {$basha[$i] = $Line[3]}
    if ($Line[0] eq "Background" && $Line[2] eq "(averin):") {$baver[$i] = $Line[3]}
    if ($Line[0] eq "Averin" && $Line[1] eq "signal,") {$ares[$i] = $Line[4]}
    if ($Line[0] eq "Detector" && $Line[1] eq "pulse") {$ampd[$i] = $Line[6]}
    if ($Line[0] eq "Amplifier" && $Line[1] eq "pulse") {$ampa[$i] = $Line[6]}
    if ($Line[0] eq "Pulse" && $Line[2] eq "(deconvolved):") {$pusid[$i] = $Line[3]}
    if ($Line[0] eq "Pulse" && $Line[2] eq "(shaped):") {$pusis[$i] = $Line[3]}
    if ($Line[0] eq "Deconvolved" && $Line[2] eq "number") {
      $red[$i] = $Line[5];
      $reda[$i] = $Line[6];
    }
    if ($Line[0] eq "Deconvolved" && $Line[2] eq "FWHM") {
      $resum[$i] = $Line[8];
      $resti[$i] = $Line[9];
    }
    if ($Line[0] eq "Shaped" && $Line[2] eq "number") {$res[$i] = $Line[5]}
    if ($Line[0] eq "Moving" && $Line[5] eq "samples") {$rom[$i] = $Line[10]}
    if ($Line[0] eq "FADC" && $Line[1] eq "delay") {
      $fat[$i] = $Line[5];
      $phs[$i] = $Line[6];
      $taf[$i] = $Line[9];
      print "$i $rami[$i] $fat[$i] $rom[$i] $tisi[$i]\n";
    }
# Data for derivative normalization
    if ($Line[0] eq "Moving" && $Line[4] eq "FADC") {$fani[$i] = $Line[9]}
    if ($Line[0] eq "Derivative" && $Line[1] eq "integral:") {$deni[$i] = $Line[2]}
  }
  close(Refile);
  printf Difile "%5d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", 
#         $i,$fat[$i],$rami[$i],$rom[$i],$dam[$i],$resum[$i],$res[$i],$red[$i],$reda[$i],$resti[$i],$tisi[$i],$badec[$i],$basha[$i];
         $i,$fat[$i],$ram[$i],$phs[$i],$dam[$i],$red[$i],$res[$i],$pusid[$i],$pusis[$i],$badec[$i],$basha[$i],$fani[$i],$deni[$i],$ampd[$i],$ampa[$i],$rom[$i],$taf[$i],$tar,$sar,$dar;
  $suma = $suma + $res[$i];
  $sumt = $sumt + $tish[$i];
  if ($res[$i] > $maxram) {
    $maxram = $res[$i];
    $maxtim = $fat[$i];
  }
  if ($res[$i] < $minram) {
    $minram = $res[$i];
    $mintim = $fat[$i];
  }
  if ($tish[$i] > $maxti) {
    $maxti = $tish[$i];
    $maxtime = $fat[$i];
  }
  if ($tish[$i] < $minti) {
    $minti = $tish[$i];
    $mintime = $fat[$i];
  }
}
close (Difile);

print "Minimal signal amplitude: $minram $mintim\n";
print "Maximal signal amplitude: $maxram $maxtim\n";

print "Minimal signal arrival time and its phase: $minti $mintime\n";
print "Maximal signal arrival time and its phase: $maxti $maxtime\n";

#print "suma $suma\n";

$ramean = $suma/$ni;
$timean = $sumt/$ni;

$suma = 0;
for ($i=0;$i<$ni;$i++) {
  $suma = $suma + ($res[$i]-$ramean)*($res[$i]-$ramean);
}
#print "suma sigma $suma\n";
$simean = sqrt($suma/($ni-1));

print "Shaped signal amplitude, mean and sigma: $ramean $simean\n";

$suma = 0;
for ($i=0;$i<$ni;$i++) {
  $suma = $suma + ($tish[$i]-$timean)*($tish[$i]-$timean);
}
#print "suma sigma $suma\n";
$simean = sqrt($suma/($ni-1));

print "Shaped signal arrival time, mean and sigma: $timean $simean\n";
