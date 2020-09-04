#! /usr/bin/perl
#  Generate cachefiles and jobfiles for segments between t1 and t2.
#  There is only one job.

$ifo1 = "H";
$ifo2 = "L";

$t1 = 793177024;
$t2 = 793178048;

$type1="RDS_R_L3";  #LHO
$type2="RDS_R_L3";  #LLO

$list1=`ligo_data_find -o $ifo1 -t $type1 -s $t1 -e $t2`;
$list1 =~ s/gsiftp.*\n//g;
$list1 =~ s/.*data.*\n//g;
$list1 =~ s/file:\/\/localhost//g;
$gps1=$list1;
$gps1 =~ s/.*$type1-(.*)-[0-9]*\.gwf\n/\1\n/g;
if ($list1 =~ m/No files found.*/) {print "No files found for $ifo1";}

$list2=`ligo_data_find -o $ifo2 -t $type2 -s $t1 -e $t2`;
$list2 =~ s/gsiftp.*\n//g;
$list2 =~ s/.*data.*\n//g;
$list2 =~ s/file:\/\/localhost//g;
$gps2=$list2;
$gps2 =~ s/.*$type2-(.*)-[0-9]*\.gwf\n/\1\n/g;
if ($list2 =~ m/No files found.*/) {print "No files found for $ifo2";}

open(gpsfile1,">cachefiles/gpsTimes$ifo1.1.txt");
open(gpsfile2,">cachefiles/gpsTimes$ifo2.1.txt");
open(framefile1,">cachefiles/frameFiles$ifo1.1.txt");
open(framefile2,">cachefiles/frameFiles$ifo2.1.txt");
print framefile1 $list1;
print framefile2 $list2;
print gpsfile1 $gps1;
print gpsfile2 $gps2;

close(gpsfile1); close(gpsfile2); close(framefile1); close(framefile2);

open(jobfile,">jobfiles/test_jobs.txt");
$diff=$t2 - $t1;
print jobfile "%seg    start      stop    duration\n";
print jobfile "   1  $t1  $t2   $diff\n";
