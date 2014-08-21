#!/usr/bin/perl
 
$sdir = $ARGV[0];

$cgrid = "xrun_";
$cpbs = ".pbs";
$cl = "/";

$gridfile = $cgrid.$sdir.$cpbs;
open(Pbsfile,'>',$gridfile);
  printf Pbsfile "#!/bin/tcsh\n";
  printf Pbsfile "#PBS -l walltime=168:00:00,nodes=1:core2duo,mem=1800MB,vmem=1800MB\n";
  printf Pbsfile "#PBS -q pub7d\n";
  printf Pbsfile "cd /work/1/lubinski/sy/spectra/$sdir\n";
  printf Pbsfile "source /app/lheasoft/setup-app-lheasoft-6.11.csh\n";
  printf Pbsfile "xspec - grid.xcm\n";
close(Pbsfile);
system "qsub ./$gridfile";  
