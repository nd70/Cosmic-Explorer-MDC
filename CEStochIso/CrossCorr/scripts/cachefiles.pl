#! /usr/bin/perl
# E. Thrane on June 4, 2010
# Creates cachefiles from $jobfile.
# Edit the part of the script enclosed by dashed lines to create a different
# set of cachefiles

# S5H1L1 parameters------------------------------------------------------------
$jobfile = "S5H1L1_full_run.txt";
$type1 = "RDS_R_L1";
$type2 = "RDS_R_L1";
$IFO1 = "H1";
$IFO2 = "L1";
# -----------------------------------------------------------------------------

# Create cachfiles
open(jobfile,$jobfile);
$site1 = $IFO1;
$site2 = $IFO2;
$site1 =~ s/(.)./\1/;
$site2 =~ s/(.)./\1/;
# job number
$i = 1;

# make cachefile/ dir if necessary
if (!(-d "./cachefiles/")) {system "mkdir cachefiles";}

while(<jobfile>) {                                      # Loop over jobs
    # Create frameFiles and gpsTime files
    $framefile1 = "cachefiles/frameFiles$site1.$i.txt";
    $framefile2 = "cachefiles/frameFiles$site2.$i.txt";
    open(framefile1,">$framefile1");  open(framefile2,">$framefile2");
    $gpsfile1 = "cachefiles/gpsTimes$site1.$i.txt";
    $gpsfile2 = "cachefiles/gpsTimes$site2.$i.txt";
    open(gpsfile1,">$gpsfile1");  open(gpsfile2,">$gpsfile2");
    # Parse jobfile and obtain frames using ligo_data_find
    # format: job# start stop duration
#    $start = $_;  $start =~ s/.* (.*) (.*) .*/\1/;  chomp($start);
#    $stop = $_;   $stop =~ s/.* (.*) (.*) .*/\2/;   chomp($stop);
    chomp($_);
    $start = $_;
    $start =~ s/ *([0-9]*) *([0-9]*) *([0-9]*) *([0-9]*) */\2/;
    $stop = $_;
    $stop =~ s/ *([0-9]*) *([0-9]*) *([0-9]*) *([0-9]*) */\3/;

# diagnostic test
#print "ligo_data_find -o $site1 -t $type1 -s $start -e $stop -u 'file'\n";
#    exit;

    $frames1=`ligo_data_find -o $site1 -t $type1 -s $start -e $stop -u 'file'`;
    $frames1 =~ s/file:\/\/localhost//g;
    $frames2=`ligo_data_find -o $site2 -t $type2 -s $start -e $stop -u 'file'`;
    $frames2 =~ s/file:\/\/localhost//g;

    if ($frames1=~m/No files found!.*/) {print "Error with $site1 $start\n";}
    if ($frames2=~m/No files found!.*/) {print "Error with $site2 $start\n";}

    # Extract GPS times
    $times1=$frames1;  $times2=$frames2;
    $times1=~s/.*$type1-(.*)-.*/\1/g;  $times2=~s/.*$type2-(.*)-.*/\1/g;
 
   # Print data to cachefiles
    print framefile1 "$frames1";  print framefile2 "$frames2";
    print gpsfile1 "$times1";  print gpsfile2 "$times2";
    # Clean up
    close(framefile1);  close(framefile2);
    close(gpsfile1);  close(gpsfile2);
    $i = $i+1;
}
close(jobfile);
