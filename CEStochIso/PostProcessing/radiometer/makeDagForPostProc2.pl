#! /usr/bin/env perl
# Created by Meyers
# meyers@physics.umn.edu
#######################################
# Creates an executable and submit script
# for condor for radiometer
# postprocessing
# Arguments should be:
# 0 --------> analysis summary file name (see README.txt)
# 1 --------> post processing file name
# 2 --------> number of jobs to combine in single postproc
#             job. typically set to 250
#######################################

use strict;
use warnings;
use File::Basename;
# read in arguments
my $analysisFile;
my $post_proc_params_file;
my $jobsToCombine;
my $dagName;
my $exeDir;
my $cjobs;
######
# check first arg
######
if (defined $ARGV[0]){
  $analysisFile = $ARGV[0];
  print "Summary of param and job files is in: \n$analysisFile\n\n";
  }
else{die "Need a first argument\n";}
########
# check second arg
########
if (defined $ARGV[1]){
  $post_proc_params_file = $ARGV[1];
  unless (-e $post_proc_params_file){die "Post proc params file: $post_proc_params_file  doesn't exist!\n"}
  else{print "Post Proc Params file exists: $post_proc_params_file\n\n"}
}
else{die "Need a second argument specifying post processing script\n";}

##########
# check third arg
# if it doesn't exist
# set it to 250
##########
if (defined $ARGV[2]){
   $jobsToCombine = $ARGV[2]; 
   print "Jobs To combine: $jobsToCombine\n";
   }
else{
  $jobsToCombine = 250;
  print "Jobs To Combine not defined. Setting it to 250\n";
}

#######
# check if dag name is defined
# if not set it automatically
# to rad_pproc
#######
if (defined $ARGV[3]){$dagName = $ARGV[3];
   print "dag is named: $dagName\n";
}
else{
  $dagName = "rad_pproc";
  print "dag name is not specified.\nautomatically setting it to \"./rad_pproc\".\n";
}

$exeDir = $exeDir."./";
print "executable directory: $exeDir\n";


# read in to a 2D array the necessary information
# for ifo pair, param files and job files used
open(INPUTS,"<".$analysisFile);
my @analysisINFO; 
while (my $line = <INPUTS>) {

  push @analysisINFO,[split(' ',$line)] unless $line eq "\n";

}
close(INPUTS);

my($no, $dagDir,$no2) = fileparse($dagName);

   unless(-e $dagDir."/out"){system("mkdir -p ".$dagDir."/out");}
   unless(-e $dagDir."/err"){system("mkdir -p ".$dagDir."/err");}
   

################################
# first line of file should give
# location of params files first
# then jobs files for this analysis
#################################

# how many are in the analysis file?
my $num = scalar(@analysisINFO);
# for keeping track of jobs
my $m=1;

# dag file
open(DAG,">$dagName.dag");

print "dag file is: $dagName.dag\n";
####################################
# loop over rows in analysis summary
# file.
####################################
for(my $j=0;$j<$num;$j++){

   # get the name of job file
   my $file = $analysisINFO[$j][2];
###################################################
# figure out how many combined jobs I need based on
# number of stochastic jobs
###################################################
  #open job file 
  open(JOBS,"<".$file) or die "Can't open: $!";
  my @job_array = <JOBS>; 
  my $jobMax    = scalar(@job_array);
   use POSIX;
   $cjobs = ceil($jobMax/$jobsToCombine);
  close(JOBS);

# check that params file exists..
unless(-e $analysisINFO[$j][1]){die "$analysisINFO[$j][1] doesn't exist!\n";}

# write to dag file
  for(my $k=1;$k<=$cjobs;$k++){
    print DAG "JOB $m $dagName.sub\n";
    print DAG "VARS $m ";
    print DAG "paramsFile=\"".$analysisINFO[$j][1]."\" ";
    print DAG "postProcParamsFile=\"$post_proc_params_file\" ";
    print DAG "jobsFile = \"".$analysisINFO[$j][2]."\" ";
    print DAG "epoch = \"".$analysisINFO[$j][0]."\" ";
    print DAG "jobNumber=\"$k\" ";
    print DAG "num=\"".$m."\" \n\n";
    $m = $m + 1;
  }
}

#print sum local rad job
print DAG "JOB $m $dagName.sub\n";
print DAG "VARS $m ";
print DAG "paramsFile=\"placeholder\" ";
print DAG "postProcParamsFile=\"$post_proc_params_file\" ";
print DAG "jobsFile = \"placeholder\" ";
print DAG "epoch = \"SumLocalRad\" ";
print DAG "jobNumber=\"1\" ";
print DAG "num=\"".$m."\" \n\n";

print DAG "PARENT 1 CHILD ";

for (my $k=2;$k<=$m;$k++){
  print DAG "$k ";
}

print DAG "\n";
print DAG "PARENT ";

for (my $k=1;$k<$m;$k++){
  print DAG "$k ";
}
print DAG "CHILD $m\n";
  
close(DAG); # close dag file

# write submit file
open(DSUB,">$dagName.sub");
  print "sub file is: $dagName.sub\n";
  print DSUB "Executable     = ".$exeDir."radiometer_postproc\n";
  print DSUB "Output         = ".$dagDir."out/\$(num).out\n";
  print DSUB "Error          = ".$dagDir."err/\$(num).err\n";
  print DSUB "Log            = ".$dagDir."out/log.log\n";
  print DSUB "Arguments      = \$(paramsFile) \$(postProcParamsFile) \$(jobsFile) \$(epoch) \$(jobNumber)\n";
  print DSUB "request_memory = 4000\n";
  print DSUB "Universe       = vanilla\n";
  print DSUB "Notification   = never\n";
  print DSUB "environment    = LD_LIBRARY_PATH=/usr/local/MCR/matlab_r2013a/v81/runtime/glnxa64:/usr/local/MCR/matlab_r2013a/v81/bin/glnxa64:/usr/local/MCR/matlab_r2013a/v81/sys/os/glnxa64:/usr/local/MCR/matlab_r2013a/v81/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/usr/local/MCR/matlab_r2013a/v81/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/local/MCR/matlab_r2013a/v81/sys/java/jre/glnxa64/jre/lib/amd64:\${LD_LIBRARY_PATH};XAPPLRESDIR=/usr/local/MCR/matlab_r2013a/v81/X11/app-defaults\n";
  print DSUB "getenv         = true\n";
  print DSUB "+MaxHours      = 12\n";
  print DSUB "Queue 1";
close(DSUB); # close submit file
# EOF
