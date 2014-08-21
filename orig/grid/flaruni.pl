#!/usr/bin/perl
 
$sdir = $ARGV[0];

$cgrid = "flaruni_";
$cpbs = ".pbs";
$cl = "/";

$gridfile = $cgrid.$sdir.$cpbs;
open(Pbsfile,'>',$gridfile);
  printf Pbsfile "#!/bin/tcsh\n";
  printf Pbsfile "#PBS -l walltime=168:00:00,nodes=1,mem=1500MB,vmem=1500MB\n";
  printf Pbsfile "#PBS -q pub7d\n";
  printf Pbsfile "cd /work/1/lubinski/ctanew/signal/$sdir\n";
  printf Pbsfile "./flaszka.pl 100000\n";
close(Pbsfile);
system "qsub ./$gridfile";  
