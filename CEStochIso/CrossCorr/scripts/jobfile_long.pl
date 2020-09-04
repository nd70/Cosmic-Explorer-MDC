#! /usr/bin/perl
# E. Thrane on May 31, 2010 - creates mini jobfiles for long jobs
# After this, run combine_jobfiles.pl to make a single jobfile.

# User-defined variables ------------------------------------------------------
$veto_path = "/archive/home/ethrane/S6";
$jobdir = "./jobfiles";
$IFO1 = "H1";
$IFO2 = "L1";
$type1 = "H1_LDAS_C00_L2";
$type2 = "L1_LDAS_C00_L2";
# exclude jobs shorter than $min_dur
$min_dur = 156;      
$START_GPS = 959006522;
$END_GPS = 959067722;
$DeltaT = 3600*10;
#-----------------------------------------------------------------------------

$JOB = 1;
for ($startGPS=$START_GPS; $startGPS<=$END_GPS; $startGPS=$startGPS+$DeltaT) {
    # Update endGPS
    $endGPS = $startGPS+$DeltaT;
    # Define veto flags in ??_vetoflags.txt files.
    $exclude1 = `cat $veto_path/$IFO1"_vetoflags.txt"`;
    $exclude2 = `cat $veto_path/$IFO2"_vetoflags.txt"`;
    $exclude = $exclude1.",".$exclude2;
    # Create jobfile_$JOB.xml
    system "/opt/lscsoft/glue/bin/ligolw_segment_query --query-segments -d -s $startGPS -e $endGPS --include-segments '$IFO1:DMT-SCIENCE:1,$IFO2:DMT-SCIENCE:1' --exclude-segments '$exclude' > $jobdir/jobfile_$JOB.xml";
    # Parse jobfile_$JOB.xml to create jobfile_$JOB.tmp
    system "/opt/lscsoft/glue/bin/ligolw_print $jobdir/jobfile_$JOB.xml --table segment --column start_time --column end_time --delimiter \" \" > $jobdir/jobfile_$JOB.tmp";
    # Remove jobfile_$JOB.xml
    system "rm $jobdir/jobfile_$JOB.xml";
    # Add columns to jobfile_$JOB
    open(tmpfile,"$jobdir/jobfile_$JOB.tmp");
    open(jobfile,">$jobdir/jobfile_$JOB");
    $i=1;
    $tot_dur=0;
    while (<tmpfile>) {
        $start = $_;  $start =~ s/(.*) (.*)/\1/;  chomp($start);
        $stop = $_;   $stop =~ s/(.*) (.*)/\2/;   chomp($stop); 
        $dur = $stop-$start;
        if ($dur>=$min_dur) {
            $tot_dur = $tot_dur + $dur;
            print jobfile "$i $start $stop $dur\n";
            $i = $i+1;
        }
    }
    close(tmpfile);  close(jobfile);
    system "rm $jobdir/jobfile_$JOB.tmp";
    
    $JOB=$JOB+1;
}
