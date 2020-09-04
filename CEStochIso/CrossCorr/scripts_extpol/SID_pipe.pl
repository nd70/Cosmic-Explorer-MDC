#! /usr/bin/perl
#
# usage: SID_pipe.pl paramfile.txt jobfile.txt dagfile.dag
# Please include full path names. -ET

#Open paramfile and determine SID frame cache path and frame type.
$paramfile = shift;  $jobfile = shift;  $dagfile = shift;
if ( !($paramfile =~ m/.*\.txt/) && !($jobfile =~ m/.*\.txt/) 
     && !($dagfile =~ m/.*\.dag/) ) {
    print " usage: SID_pipe.pl paramfile.txt jobfile.txt\n";
    exit;
}

open(parm,$paramfile);
while (<parm>) {
    if($_ =~ m/intFrameCachePath.*/) {
	$intFrameCachePath = $_;
	$intFrameCachePath =~ s/intFrameCachePath (.*)/\1/;
	chomp($intFrameCachePath);
    }
    if($_ =~ m/frameType1.*/) {
	$frameTypeInt = $_;
	$frameTypeInt =~ s/frameType1 (.*)/\1/;
	chomp($frameTypeInt);
	$BL = $frameTypeInt;  $BL =~ s/(.*)-.*/\1/;
	$type = $frameTypeInt;  $type =~ s/.*-(.*)/\1/;
    }
}
close(parm);
if (!(-e $intFrameCachePath)) {
    print "intFrameCachePath does not exist: please modify paramfile.\n";
    exit;
}


#Open jobfile, for each line (segment), determine gpsStart/End times.
#For each segment, call LSCDataFind and write frame and gps cachefiles.
open(jobs,$jobfile);  open(dag,">$dagfile");
$n=1; #job number
#Define environmental variables.
$LD_LIBRARY_PATH = "$ENV{'LD_LIBRARY_PATH'}";  $HOME = "$ENV{'HOME'}";
$nMax=1E99;       #diagnostic tool
#$nMax=24;        #create only $nMax jobs
while(<jobs>) {
    if($n<=$nMax) { #------------------------------------------------------
	#write to dagfile
	print dag "JOB $n stochastic_pipe.sub\n";
	print dag "VARS $n paramsFile=";
	print dag '"';  print dag "$paramfile";
	print dag '" jobsFile="';  print dag "$jobfile";
	print dag '" jobnumber="';  print dag "$n";    
	print dag '" home="';  print dag "$HOME";
	print dag '" ld_library_path="';  print dag "$LD_LIBRARY_PATH";
	print dag '"';  print dag "\n\n";
	
	#write to jobfile
	$gpsStart = $_;  $gpsEnd = $_;
	$gpsStart =~ s/ *[0-9]* *([0-9]*) .*/\1/;  chomp($gpsStart);
	$gpsEnd =~ s/ *[0-9]* *[0-9]* *([0-9]*) .*/\1/;  chomp($gpsEnd);
	$frameFile = $intFrameCachePath."frameFiles".$BL.".".$n.".txt";
	$gpsTimes = $intFrameCachePath."gpsTimes".$BL.".".$n.".txt";
	if ((-e $frameFile) && (-e $gpsTimes)) {}
	else {
	    $cachefiles = `LSCdataFind -o $BL -t $type -s $gpsStart -e $gpsEnd`;
	    $cachefiles =~ s/file:\/\/localhost//g;
	    $cachefiles =~ s/gsiftp.*\n//g;
	    open(Ffile,">$frameFile");  print Ffile $cachefiles;  close(Ffile);
	    $cachefiles =~ s/.*-([0-9]*)-[0-9]*\.gwf/\1/g;
	    open(gfile,">$gpsTimes");  print gfile $cachefiles;  close(gfile);
	}
	print "$n ";  
	$n++;
    } #----------------------------------------------------------------
}

#finish up dagfile
print dag "JOB 0 stochastic_pipe.sub\n";
print dag "VARS 0 paramsFile=";
print dag '"';  print dag "$paramfile";
print dag '" jobsFile="';  print dag "$jobfile";
print dag '" jobnumber="';  print dag "0";    
print dag '" home="';  print dag "$HOME";
print dag '" ld_library_path="';  print dag "$LD_LIBRARY_PATH";
print dag '"';  print dag "\n\n";
print dag "PARENT 0 CHILD ";
for ($i=1; $i<$n; $i++) {print dag "$i ";}  print dag "\n\n";

#close files
close(jobs);  close(dag);  print "\n";
