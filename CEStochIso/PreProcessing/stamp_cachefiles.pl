#! /usr/bin/perl
# E. Thrane:  Define the variables at the head of this script.
# STAMP cachefiles will be created in cachefiles/.  If you want to do the
# entire jobfile, set $maxjobs to a very large number.
#
# location of preproc frames
$FrDir = "/archive/home/ethrane/preproc/frames";
# frame names of the form $str-$gps-$dur.gwf
#$str = "HL-noWindow";
#$str = "S5H1L1";
$str = "S5H1L1_inj";
# jobfile (same format as stochastic jobfiles)
$jobfile = "S5H1L1_full_run.txt";
# sites is a two letter code for the IFO locations.  HL=Hanford-Livingston.
$sites = "HL";
# Max number of jobs for which you will create cachefiles.
$maxjobs = 20;

# make cachefile/ dir if necessary
if (!(-d "./cachefiles/")) {system "mkdir cachefiles";}

$i = 1;
open(jobfile, $jobfile);
while(<jobfile>) {                      # Loop over jobs
    if ($i>$maxjobs) {exit;}

    # Create frameFiles and gpsTime files
    $framefile = "cachefiles/frameFiles$sites.$i.txt";
    open(framefile,">$framefile");
    $gpsfile = "cachefiles/gpsTimes$sites.$i.txt";
    open(gpsfile,">$gpsfile");

    # Parse jobfile
    $start = $_;  $stop = $_;
    $start =~ s/ *[0-9]* *([0-9]*) *([0-9]*) *[0-9]* */\1/;
    $stop =~ s/ *[0-9]* *([0-9]*) *([0-9]*) *[0-9]* */\2/;
    chomp($start);  chomp($stop);
    print "$start - $stop\n";

    foreach $file (<$FrDir/$str*>) {
	$gps = $file;
	$gps =~ s/$FrDir\/$str-(.*)-[0-9]*.gwf/\1/;
	print "$gps\n"; 	print "$file\n";
	if ($gps>=$start && $gps<=$stop) {
          # Print data to cachefiles
	  print framefile "$file\n";
	  print gpsfile "$gps\n";  
	}
    }

    # Clean up
    close(framefile);  close(gpsfile);
    $i = $i+1;
}
close(jobfile);
