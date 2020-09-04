#! /usr/bin/env perl
# T. Prestegard
# Creates executable and condor submit script.
# All jobs go in the same cluster on condor.

use strict;
use warnings;

#############################################################
# Parameters
#############################################################
# Set parameters.
#my $home = $ENV{'HOME'};
#my $preproc = $home."/matapps/packages/stochastic/trunk/PreProcessing";
my $home = "/home/sgwynne.crowder";

# Parameters to change.
#my $name = "preproc"; # Name of this directory in ~/condor.
#my $paramfile = $preproc."/Xi_bs/params.txt";
my $paramfile = "/home/sgwynne.crowder/sgwb/O1_zerolag/C02/input/paramfiles/C02_params_a0_fullComb.txt";
#my $jobfile = $preproc."/Xi_bs/new_padded_jobfile_1000.txt";
my $jobfile = "/home/sgwynne.crowder/sgwb/O1_zerolag/C02/input/jobfiles/JOB-FILE-1126623617-1136649617-master-C02-lowf.dat";
#my $exe = "preproc"; # Name of executable (compiled Matlab code).
my $n1 = 1; # First job in the jobfile that you want to run.
my $nMax = 2000; # Set this to a very large number if you want all jobs in the jobfile.

# Parameters that do not need to be changed.
#my $dir = $home."/condor/".$name;
#my $cexec = $name."_exec.sh";

#############################################################
# Code.
#############################################################

# Get number of jobs in jobfile.
open(JOBS,"<".$jobfile);
my @job_array = <JOBS>; close(JOBS);
my $jobMax = (scalar(@job_array) >= $nMax) ? $nMax : scalar(@job_array);
close(JOBS);

# Make executable file.
#open(SH,">".$dir."/sh/".$cexec);
#print SH "#! /usr/bin/env bash\n\n";
#print SH "cd ".$home."/sgwb/trunk/S5/matlab\n";
#print SH "source matlab_script_64bit.sh\n";
#print SH "cd ".$preproc."\n";
#print SH "./".$exe." \$1 \$2 \$3\n\n";
#print SH "# EOF\n";
#close(SH);
#qx(chmod u+x $dir/sh/$cexec);

# Make dagfile
#open(DAG,">".$dir."/sh/".$name.".dag");
open(dag,">./a0_fullComb.dag");

my $j=1;
for (my $i=$n1; $i<=$jobMax; $i=$i+1) {
#    print DAG "JOB ".$j." ".$name.".sub\n";
    print dag "JOB $j /home/sgwynne.crowder/sgwb/O1_zerolag/C02/input/command/a0_fullComb.sub\n";
#    print DAG "VARS ".$j." paramfile=\"".$paramfile."\" jobfile=\"".$jobfile."\" jobNumber=\"".$i."\"\n\n";
     print dag "VARS ".$j." paramsFile=\"".$paramfile."\" jobsfile=\"".$jobfile."\" jobNumber=\"".$i."\"\n\n";
#    print dag "VARS $j ";
#    print dag "paramsFile=\"/home/sgwynne.crowder/sgwb/O1_timeshift/week1/input/paramfiles/H1L1_O1_150918_150926_params_deltasigmatest_20-1639Hz_a5.txt\" ";
#    print dag "jobsfile=\"/home/sgwynne.crowder/sgwb/O1_timeshift/week1/input/jobfiles/JOB-FILE-1126623617-1127271617.dat\" ";
#    print dag "jobNumber=\"$j\" ";
#    print dag "\n\n";
    $j = $j + 1;
}

print dag "PARENT 1 CHILD ";
for (my $i=2; $i<$j; $i=$i+1) { print dag $i." "; }
print dag "\n";
close(dag);

# Make dag submit file.
#open(DSUB,">".$dir."/sh/".$name.".sub");
#print DSUB "Executable   = ".$dir."/sh/".$cexec."\n";
#print DSUB "Output       = ".$dir."/out/\$(jobNumber).out\n";
#print DSUB "Error        = ".$dir."/err/\$(jobNumber).err\n";
#print DSUB "Log          = /usr1/prestegard/".$name."_dag.log\n";
#print DSUB "Arguments    = \$(paramfile) \$(jobfile) \$(jobNumber)\n";
#print DSUB "Requirements = Memory >= 128 && OpSys == \"LINUX\"\n";
#print DSUB "Universe     = vanilla\n";
#print DSUB "Notification = never\n";
#print DSUB "getenv       = True\n";
#print DSUB "+MaxHours    = 12\n";
#print DSUB "Queue 1";
#close(DSUB);

# EOF
