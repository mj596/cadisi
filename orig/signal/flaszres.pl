#!/usr/bin/perl

# FlasCam simulation tools

use POSIX;

$infile = $ARGV[0];
$ico = $ARGV[1];

open (Ifile,'<',$infile);

$i = 0;
$minval = 1000000;
$maxval = -100000;
$suma = 0;
while (<Ifile>) {
  chomp;
  @Line = split;
  if ($Line[1] < 0.0001 && $Line[2] < 0.0001) {next}
  $val[$i] = $Line[$ico];
  $suma = $suma + $val[$i];
  if ($val[$i] > $maxval) {
    $maxval = $val[$i];
  }
  if ($val[$i] < $minval) {
    $minval = $val[$i];
  }
  $i++;
}
close (Ifile);
$ni = $i;

print "$ni lines, minimal and maximal value: $minval $maxval\n";

$vamean = $suma/$ni;

$suma = 0;
for ($i=0;$i<$ni;$i++) {
  $suma = $suma + ($val[$i]-$vamean)*($val[$i]-$vamean);
}
#print "suma sigma $suma\n";
$simean = sqrt($suma/($ni-1));

print "Mean and sigma: $vamean $simean\n";
