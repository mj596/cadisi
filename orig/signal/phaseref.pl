#!/usr/bin/perl

# Find phase for which both deconvolved and shapes signal are closest to 10

use POSIX;

$infile = $ARGV[0];

open(Ifile,'<',$infile);

$dem = 1.0;
$det = 1.0;
while (<Ifile>) {
  chomp;
  @Line = split;
  $rd = $Line[5];
  $rs = $Line[6];
  $td = $Line[7];
  $ts = $Line[8];
  $de = abs($rs-10)+abs($rd-10);
  $dt = abs($td-905.444)+abs($ts-908.2);
  if ($de < $dem) {
    printf "$de $Line[1] $Line[3] $rd $rs\n";
    $dem = $de;
  }
  if ($dt < $det) {
    printf "$dt $Line[1] $Line[3] $td $ts\n";
    $det = $dt;
  }
}
close(Ifile);
