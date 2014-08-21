#!/usr/bin/perl

# FlasCam simulation tools
# Compute mean and sigma for a sum of two distributions, e.g. 10 and 20 p.e. for 30 p.e. test

use POSIX;

$infile1 = $ARGV[0];
$infile2 = $ARGV[1];
$ico = $ARGV[2];

open (Ifile,'<',$infile1);

$i = 0;
while (<Ifile>) {
  chomp;
  @Line = split;
  if ($Line[1] < 0.0001 && $Line[2] < 0.0001) {next}
  $val1[$i] = $Line[$ico];
  $i++;
}
close (Ifile);
$n1 = $i;

open (Ifile,'<',$infile2);

$i = 0;
while (<Ifile>) {
  chomp;
  @Line = split;
  if ($Line[1] < 0.0001 && $Line[2] < 0.0001) {next}
  $val2[$i] = $Line[$ico];
  $i++;
}
close (Ifile);
$n2 = $i;

print "$n1 $n2 lines\n";

if ($n1 >= $n2) {$nn = $n2}
if ($n1 < $n2) {$nn = $n1}

$minval = 1000000;
$maxval = -100000;
$suma = 0;
for ($i=0;$i<$nn;$i++) {
  $val[$i] = $val1[$i]+$val2[$i];
  $suma = $suma + $val[$i];
  if ($val[$i] > $maxval) {
    $maxval = $val[$i];
  }
  if ($val[$i] < $minval) {
    $minval = $val[$i];
  }
}

print "Minimal and maximal value: $minval $maxval\n";

$vamean = $suma/$nn;

$suma = 0;
for ($i=0;$i<$nn;$i++) {
  $suma = $suma + ($val[$i]-$vamean)*($val[$i]-$vamean);
}
#print "suma sigma $suma\n";
$simean = sqrt($suma/($nn-1));

print "Mean and sigma: $vamean $simean\n";
