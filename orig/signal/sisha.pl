#!/usr/bin/perl

# Signal shape transformations
# Quoted pages of the FlashCam report, version 1.9, 03.04.2012

use POSIX;

# Time interval, ns
$x1 = 0;
$x2 = 200;

# Pulse arrival time, ns
$x0 = 100;

# Time resolution
$dx = 0.1;

# Time bins
$nx = int(($x2-$x1)/$dx+0.5);

# Convolution interval
$c1 = -100;
$c2 = 100;

# Convolution bins
$nc = int(($c2-$c1)/$dx+0.5);

# Original PMT pulse shape, page 36
# Rise time
$tir = 1.5;
# Fall time adjusted to get FWHM = 3.4
$tif = 2.05;
# Amplitude
$ao = 910;

open (Osfile,'>','oshape.dat');
open (Osdfile,'>','odshape.dat');
open (Osffile,'>','ofshape.dat');
open (Osifile,'>','oishape.dat');
open (Scofile,'>','conshape.dat');
open (Cdefile,'>','codeshape.dat');
open (Avefile,'>','aveshape.dat');
open (Davfile,'>','davshape.dat');
open (Dfifile,'>','dfishape.dat');
open (Decfile,'>','decoshape.dat');

# Hybrid gaussian instead
# Parameters fitted with fithygauss
# 250 MHz signal
$ao = 284.18;
$xpa = $x0 + 2.95;
$patif = 1.4154;
$tau = 0.83725;
# 1 GHz signal
$ao = 143.12;
$xpa = $x0 + 3.01;
$patif = 1.425;
$tau = 0.90954;
# 1 GHz input
#$ao = 100.79;
#$xpa = $x0 + 4.94;
#$patif = 2.0403;
#$tau = 0.90141;

$osmax = 0;
for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  if ($x[$i] < $x0) {$u = 0}
  if ($x[$i] >= $x0) {$u = 1}
#  $ex1 = exp(-($x[$i]-$x0)/$tir);
#  $ex2 = exp(-($x[$i]-$x0)/$tif);
#  $os[$i] = $ao*(1-$ex1)*$ex2*$u;
  $fas = ($x[$i]-$xpa)*$tau;
  $ex2 = 0;
  if (2*$patif*$patif + $fas > 0) {
    $ex2 = exp(-($x[$i]-$xpa)*($x[$i]-$xpa)/(2*$patif*$patif + $fas));
  }
#  $pare[$i] = $apa*(1-$ex1)*$ex2*$u;
  $os[$i] = $ao*$ex2;
  if ($os[$i] > $osmax) {
    $osmax = $os[$i];
    $oxmax = $x[$i];
  }
#  print "$i $x[$i] $os[$i] $u\n";
  printf Osfile "%.1f %.4g\n", $x[$i],$os[$i];
}
close (Osfile);

print "Original pulse maximum: $oxmax $osmax\n";

$xhalf1 = 0;
$xhalf2 = 0;
$ossum = 0;
for ($i=1;$i<$nx;$i++) {
  if ($os[$i-1] < 0.5*$osmax && $os[$i] >= 0.5*$osmax) {
    $xhalf1 = $x[$i] - 0.5*($x[$i]-$x[$i-1]);
  }
  if ($os[$i-1] >= 0.5*$osmax && $os[$i] < 0.5*$osmax) {
    $xhalf2 = $x[$i] - 0.5*($x[$i]-$x[$i-1]);
  }
  $ossum = $ossum + $dx*$os[$i];
}
$fwhm = $xhalf2-$xhalf1;

print "Original pulse FWHM: $fwhm $xhalf1 $xhalf2\n";
print "Original pulse integral: $ossum\n";

# Original pulse integral

$osimax = 0;
$osisum = 0;
for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  if ($x[$i] < $x0) {$u = 0}
  if ($x[$i] >= $x0) {$u = 1}
  $cex1 = exp($x0/$tir);
  $cex2 = exp($x0/$tif);
  $ctau = $tir*$tif/($tir+$tif);
  $tex1 = exp(-$x[$i]/$tir);
  $tex2 = exp(-$x[$i]/$tif);
  $osi[$i] = $ao*$u*$cex2*($ctau*$cex1*$tex1*$tex2 - $tif*$tex2);
#  print "$x[$i] $osi[$i]\n";
  if ($x[$i] >= $x0) {$osisum = $osisum + $dx*$osi[$i]}
  if ($osi[$i] > $osimax) {
    $osimax = $osi[$i];
    $oximax = $x[$i];
  }
}

for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  $osi[$i] = $ossum*$osi[$i]/$osisum;
#  print "$i $x[$i] $os[$i] $u\n";
  printf Osifile "%.1f %.4g\n", $x[$i],$osi[$i];
}
close (Osifile);

print "Integral of original pulse: $oximax $osimax $osisum\n";

# Original pulse derivative

$osdmax = 0;
for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  if ($x[$i] < $x0) {$u = 0}
  if ($x[$i] >= $x0) {$u = 1}
  $ex1 = exp(-($x[$i]-$x0)/$tir);
  $ex2 = exp(-($x[$i]-$x0)/$tif);
  $osd[$i] = $ao*$u*$ex2*($ex1/$tir - (1-$ex1)/$tif);
  if ($osd[$i] > $osdmax) {
    $osdmax = $os[$i];
    $oxdmax = $x[$i];
  }
#  print "$i $x[$i] $os[$i] $u\n";
  printf Osdfile "%.1f %.4g\n", $x[$i],$osd[$i];
}
close (Osdfile);

print "Derivative of original pulse: $oxdmax $osdmax\n";

# Filtered original pulse: initial + derivative

# filterin constant
$afi = 8;

$osfmax = 0;
$osfsum - 0;
for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  $osf[$i] = $os[$i] - $afi*$osd[$i];
  if ($osf[$i] < 0) {$osf[$i] = 0}
  $osfsum = $osfsum + $dx*$osf[$i];
  if ($osf[$i] > $osfmax) {
    $osfmax = $osf[$i];
    $oxfmax = $x[$i];
  }
}

$osfsum1 = 0;
for ($i=0;$i<$nx;$i++) {
  $x[$i] = $x1 + $dx*$i;
  $osf[$i] = $ossum*$osf[$i]/$osfsum;
  $osfsum1 = $osfsum1 + $dx*$osf[$i];
#  print "$i $x[$i] $os[$i] $u\n";
  printf Osffile "%.1f %.4g\n", $x[$i],$osf[$i];
}
close (Osffile);

print "Filtered (derivative) pulse maximum: $oxfmax $osfmax $osfsum1\n";

# Hybrid gaussian model of averaged pulse shape
$hno = 100;
$hce = $x0 + 7.5; 
$hdt = 5.1;
$htu = 3.5;
# Fithygauss, 2-33
$hno = 79.121;
$hce = $x0 + 7.17;
$hdt = 5.2106;
$htu = 2.9575;

$havesum = 0;
for ($i=0;$i<$nx;$i++) {
  $fas = ($x[$i]-$hce)*$htu;
  $ex2 = 0;
  if (2*$hdt*$hdt + $fas > 0) {
    $hexe = exp(-($x[$i]-$hce)*($x[$i]-$hce)/(2*$hdt*$hdt + $fas));
  }
  $have[$i] = $hno*$hexe;
  $havesum = $havesum + $dx*$have[$i];
}

print "Average signal model integral: $havesum\n";

for ($i=0;$i<$nx;$i++) {
  $have[$i] = $ossum*$have[$i]/$havesum;
  printf Avefile "%.1f %.4g\n", $x[$i],$have[$i];
}
close (Avefile);

# Derivative of the hybrid gaussian model of averaged pulse shape

$havedsum = 0;
for ($i=0;$i<$nx;$i++) {
  $fas = ($x[$i]-$hce)*$htu;
  $deh = ($x[$i]-$hce);
  $hde = 2*$hdt*$hdt+$fas;
  $ex2 = 0;
  if (2*$hdt*$hdt + $fas > 0) {
    $hexe = exp(-($x[$i]-$hce)*($x[$i]-$hce)/(2*$hdt*$hdt + $fas));
    $dexe = ($htu*$deh**2-2*$deh*$hde)/($hde*$hde);
  }
  $haved[$i] = $hno*$hexe*$dexe;
  $havedsum = $havedsum + $dx*abs($haved[$i]);
}

print "Average derivative signal model integral: $havesum\n";

for ($i=0;$i<$nx;$i++) {
  $haved[$i] = $ossum*$haved[$i]/$havedsum;
  printf Davfile "%.1f %.4g\n", $x[$i],$haved[$i];
}
close (Davfile);

# Filtered signal: averaged + averaged derivative
$nave = 1.5;
$nade = 1.0;
$hafisum = 0;
for ($i=0;$i<$nx;$i++) {
  $hafi[$i] = $nave*$have[$i] + $nade*$haved[$i];
  $hafisum = $hafisum + $dx*abs($hafi[$i]);
}

print "Average filtered signal model integral: $hafisum\n";

for ($i=0;$i<$nx;$i++) {
  $hafi[$i] = $ossum*$hafi[$i]/$hafisum;
  printf Dfifile "%.1f %.4g\n", $x[$i],$hafi[$i];
}
close (Dfifile);


# Hybrid gaussian model of deconvolved pulse shape
$hno = 100;
$hce = $x0 + 4.3; 
$hdt = 3.1;
$htu = 0.02;

$hdecsum = 0;
for ($i=0;$i<$nx;$i++) {
  $fas = ($x[$i]-$hce)*$htu;
  $ex2 = 0;
  if (2*$hdt*$hdt + $fas > 0) {
    $hexe = exp(-($x[$i]-$hce)*($x[$i]-$hce)/(2*$hdt*$hdt + $fas));
  }
  $hdec[$i] = $hno*$hexe;
  $hdecsum = $hdecsum + $dx*$hdec[$i];
}

print "Deconvolved signal model integral: $havesum\n";

for ($i=0;$i<$nx;$i++) {
  $hdec[$i] = $ossum*$hdec[$i]/$hdecsum;
  printf Decfile "%.1f %.4g\n", $x[$i],$hdec[$i];
}
close (Decfile);

# Convolution and deconvolution (convolution with an inverse response)
# Simple filter response model
# Deconvolution to average , 250 MHz
$te = 4.5;
$aa = 6.5;
$bb = 1.1;
$cc = 0.0;
# Average to deconvolved, 250 MHz
$te = 1.0;
$aa = 2.0;
$bb = 2.5;
$cc = 0.0;
$dd = 0.0;
$pdesum = 0;
for ($i=0;$i<400;$i++) {
  $xu = $i*$dx;
  $yde = ($aa+$bb*$xu+$cc*$xu**2+$dd*$xu**3)*exp(-$xu/$te);
  $zde = ($bb+$cc*$xu - ($aa+$bb*$xu+$cc*$xu**2)/$te)*exp(-$xu/$te);
  $ude = $yde + $zde; 
  printf Cdefile "$xu $yde $zde $ude\n";
}
close(Cdefile);
for ($i=0;$i<$nx;$i++) {
  $oco = 0;
# t0 starts when os becomes non-zero
  if ($i > 1) {
    for ($k=0;$k<$i+1;$k++) {
      $xu = $x[$k];
#      $xu = $x[$k]-95;
      $ere1 = ($aa+$bb*$xu+$cc*$xu**2+$dd*$xu**3)*exp(-$xu/$te);
      $ere2 = ($bb+$cc*$xu - ($aa+$bb*$xu+$cc*$xu**2)/$te)*exp(-$xu/$te);
      $ere = $ere1 + $ere2; 
      $doco = $dx*$ere*$have[$i-$k];
      $oco = $oco + $dx*$ere1*$have[$i-$k+35];
    }
  }
  $osde[$i] = $oco;
#  print "$i $x[$i] $osco[$i] $oco\n";
  $pdesum = $pdesum + $osde[$i]*$dx;
}

print "Deconvolved signal integral: $pcosum\n"; 

open (Pufile,'>','pdecoshape.dat');

$pdensum = 0;
for ($i=0;$i<$nx;$i++) {
  $osde[$i] = $ossum*$osde[$i]/$pdesum;
  $pdensum = $pdensum + $osde[$i]*$dx;
  printf Pufile "%.1f %.4g\n", $x[$i],$osde[$i];
}
close (Pufile);

print "Normalized deconvolved signal integral: $pdensum\n";

# Signal filtering
open (Fofile,'>','foshape.dat');
$ftau = 0.3;
$fossum = 0;
for ($i=10;$i<$nx-1;$i++) {
  $fos[$i] = $os[$i] - $ftau*$os[$i+10];
  if ( $fos[$i] < 0 ) {$fos[$i] = 0}
  $fossum = $fossum + $dx*$fos[$i];
}
for ($i=1;$i<$nx-1;$i++) {
  $fos[$i] = $ossum*$fos[$i]/$fossum;
  printf Fofile "%.1f %.4g\n", $x[$i],$fos[$i];
}
close (Fofile);

print "Filtered pulse integral: $fossum\n";

# Preamplifier response, page 37
# Rise time
# $patir = 14;
$patir = 14;
# Fall time adjusted to get FWHM = 3.4
 $patif = 3.5;
#$patif = 2.0;
# Amplitude
 $apa = 1600;
#$apa = 100;
# Delay
# $xpa = $x0 + 6.2;
$xpa = $x0 + 6;
# Exponential gaussian hybrid
# $tau = 4.1;
$tau = 5.0;

open (Prfile,'>','prshape.dat');

$prsum = 0;
for ($i=0;$i<$nx;$i++) {
  $fas = ($x[$i]-$xpa)*$tau;
  $ex2 = 0;
  if (2*$patif*$patif + $fas > 0) {
    $ex2 = exp(-($x[$i]-$xpa)*($x[$i]-$xpa)/(2*$patif*$patif + $fas));
  }
  $pare[$i] = $apa*$ex2;
  $prsum = $prsum + $dx*$pare[$i];
}

print "Preamplifier response integral: $prsum\n";

for ($i=0;$i<$nx;$i++) {
  $pare[$i] = $ossum*$pare[$i]/$prsum;
  printf Prfile "%.1f %.4g\n", $x[$i],$pare[$i];
}
close (Prfile);

# Original pulse and preamplifier response convolution
$pasum = 0;
for ($i=0;$i<$nx;$i++) {
  $pa[$i] = $os[$i]*$pare[$i];
  $pasum = $pasum + $pa[$i]*$dx;
}

print "Preamplifier signal integral: $pasum\n"; 

open (Pafile,'>','pashape.dat');

$pansum = 0;
for ($i=0;$i<$nx;$i++) {
  $pa[$i] = $ossum*$pa[$i]/$pasum;
  $pansum = $pansum + $pa[$i]*$dx;
  printf Pafile "%.1f %.4g\n", $x[$i],$pa[$i];
}
close (Pafile);

print "Normalized preamplifier signal integral: $pansum\n";

$wex = 1.3;

# Product of two distributions
$prosum = 0;
for ($i=0;$i<$nx;$i++) {
#  $pro[$i] = $os[$i]*$pare[$i];
  $wex = 1.3 - ($x[$i]-$x0)*0.001;
  $pro[$i] = $os[$i]*1E-11*exp(($x[$i]-$x0)/$wex);
  $prosum = $prosum + $pro[$i]*$dx;
}

print "Product signal integral: $prosum\n"; 

open (Prfile,'>','proshape.dat');

$pronsum = 0;
for ($i=0;$i<$nx;$i++) {
  $pro[$i] = $ossum*$pro[$i]/$prosum;
  $pronsum = $pronsum + $pro[$i]*$dx;
  printf Prfile "%.1f %.4g\n", $x[$i],$pro[$i];
}
close (Prfile);

print "Normalized product signal integral: $pronsum\n";

# Convolution of two distributions
# Simple filter response model
# signal to input, 250 MHz
$taure = 4.5;
$aa = 6.5;
$bb = 1.2;
$cc = 0;
$dd = 0;
# signal to input, 1 GHz
$taure = 1.15;
$aa = 3.0;
$bb = 10.0;
$cc = 0;
$dd = 0;
$pcosum = 0;
for ($i=0;$i<$nx;$i++) {
  $oco = 0;
  if ($i > 1) {
    for ($k=0;$k<$i+1;$k++) {
      $xu = $x[$k];
      $etare = ($aa+$bb*$xu+$cc*$xu**2)*exp(-$xu/$taure);
      if ($i == 1700) {
        printf Scofile "$xu $etare\n";
      }
      $oco = $oco + $dx*$etare*$os[$i-$k];
    }
  }
  $osco[$i] = $oco;
  $pcosum = $pcosum + $osco[$i]*$dx;
}
close(Scofile);

print "Convolved signal integral: $pcosum\n"; 

open (Pvfile,'>','pcoshape.dat');

$pconsum = 0;
for ($i=0;$i<$nx;$i++) {
  $osco[$i] = $ossum*$osco[$i]/$pcosum;
  $pconsum = $pconsum + $osco[$i]*$dx;
  printf Pvfile "%.1f %.4g\n", $x[$i],$osco[$i];
}
close (Pvfile);

print "Normalized convolved signal integral: $pconsum\n";

# Deconvolution : convolution with an inverse response
# Simple filter response model
$taure = 3.3;
$aa = 2.5;
$bb = 2;
$cc = 0.0;
$dd = 1;
$pdesum = 0;
$etare = (1+$taure)*exp(-$taure);
$ita = int($taure/$dx);
#print "etare ita: $etare $ita\n";
for ($i=0;$i<400;$i++) {
  $xde = $i*$dx;
  $yde = ($aa+$bb*$xde+$cc*$xde**2)*exp(-$xde/$taure);
  $zde = ($bb+$cc*$xde - ($aa+$bb*$xde+$cc*$xde**2)/$taure)*exp(-$xde/$taure); 
#  printf Sdefile "$xde $zde\n";
}
#close(Sdefile);
for ($i=0;$i<$nx;$i++) {
  $oco = 0;
# t0 starts when os becomes non-zero
###  if ($i > 984) {
###    for ($k=985;$k<$i+1;$k++) {
  if ($i > 1) {
    for ($k=1;$k<$i+1;$k++) {
#      $etare = (1+($x[$k]-$x[0])/$taure)*exp(-($x[$k]-$x[0])/$taure);
#      $etare = (1+$x[$k]/$taure)*exp(-$x[$k]/$taure);
###      $xu = $x[$k]-98.5;
      $xu = $x[$k];
#      $etare = ($aa+$bb*$xu+$cc*$xu**2)*exp(-$xu/$taure);
      $etare = ($bb+$cc*$xu - ($aa+$bb*$xu+$cc*$xu**2)/$taure)*exp(-$xu/$taure); 
#      $etare = ($aa+$bb*$xu+$cc*$xu**2+$dd*$xu**3)*exp(-$xu);
#      $doco = $dx*$etare*$os[$i-$k+985];
###      $oco = $oco + $dx*$etare*$os[$i-$k+985];
      $oco = $oco + $dx*$etare*$os[$i-$k];
#      $oco = $oco + $dx*$osd[$i]*$osco[$i-$k+985];
#    if ($i == 1040) {print "$k $x[$k] $xu $etare $os[$i-$k+985] $doco $oco\n"}
    }
  }
##  $osco[$i] = $os[$i] - $oco;
  $osde[$i] = $oco;
#  print "$i $x[$i] $osco[$i] $oco\n";
  $pdesum = $pdesum + $osde[$i]*$dx;
}

print "Deconvolved signal integral: $pcosum\n"; 

#open (Pufile,'>','pdecoshape.dat');

$pdensum = 0;
for ($i=0;$i<$nx;$i++) {
  $osde[$i] = $ossum*$osde[$i]/$pdesum;
  $pdensum = $pdensum + $osde[$i]*$dx;
#  printf Pufile "%.1f %.4g\n", $x[$i],$osde[$i];
}
#close (Pufile);

print "Normalized deconvolved signal integral: $pdensum\n";

open (Infile,'<','fp37_input.dat');
$in = 0;
while (<Infile>) {
  chomp;
  @Line = split;
  $xin[$in] = $Line[0];
  $yin[$in] = $Line[1];
  $in++;
}

# Ratio of two distributions
$rasum = 0;
for ($i=0;$i<$nx;$i++) {
  if ($os[$i] > 0) {
    $yint = 0;
    for ($j=0;$j<$in-1;$j++) {
      if ($x[$i] >= $xin[$j] && $x[$i] < $xin[$j+1]) {
        $yint = $yin[$j] + ($x[$i]-$xin[$j])*($yin[$j+1]-$yin[$j])/($xin[$j+1]-$xin[$j]);
#        print "$x[$i] $os[$i] $yint\n";
      }
    }
    $pra[$i] = $yint/$os[$i];
  #  print "$x[$i] $os[$i] $yint $pra[$i]\n";
    $rasum = $rasum + $pra[$i]*$dx;
  }
}

print "Ratio signal integral: $rasum\n"; 

open (Prafile,'>','prashape.dat');

$ransum = 0;
for ($i=0;$i<$nx;$i++) {
  $pra[$i] = $ossum*$pra[$i]/$rasum;
  $ransum = $ransum + $pra[$i]*$dx;
  printf Prafile "%.1f %.4g\n", $x[$i],$pra[$i];
}
close (Prafile);

print "Normalized ratio signal integral: $ransum\n";
