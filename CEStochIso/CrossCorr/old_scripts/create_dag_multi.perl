#!/usr/bin/perl -w
# creates dag file for multi ifo analysis, using stochastic_multi.sub, 
# stochastic_multi_postproc.sub and stochastic_multi_combine.sub

use Getopt::Long;
Getopt::Long::config('no_ignore_case');

GetOptions('excluded|x:s' => \$excluded,				#x
	   'runid|i=s' => \$runid,					#i
	   'shellcomments|C' => \$shellcomments,			#C
	   'dsigmacut|K=f' => \$dsigmacut,				#K
	   'largesigmacut|L=f' => \$largesigmacut,			#L
	   'dorenormalize|z' => \$dorenormalize,			#z
	   'modifyfilter|M=s' => \$modifyfilter);				#M

# Error catching and creation of pairs (excluding any if necessary)

unless (defined($dsigmacut)){
	$dsigmacut=0.2;
}

unless (defined($largesigmacut)){
	$largesigmacut=1e20;
}

unless ($dorenormalize){
	$dorenormalize=0;
}
#if (defined($modifyfilter)){
#	unless ($dorenormalize){
#		die "Modify filter specified, but do renormalize not specified";
#	}else{
#		@filter=split(/:/,$modifyfilter);
#	}
#
#}else{
#	if ($dorenormalize){
#		die "please specify a filter";
#	}else{
#		for($i=0;$i<$npairs;$i++){
#			push(0,@filter);
#	}
#}

unless (defined($runid)) {
	die "Please give a unique identifier for your analysis";
}


if (defined($excluded)) {
	@excludedarray=split(/,/,$excluded);
	$nexcl=@excludedarray;
}else{
	$nexcl=0;
}

if (defined(@ARGV)) {
#	@ifoarray=split(/:/,$ifos);
	$nifo=@ARGV;
	$ij=0;
	for ($i = 0; $i < $nifo-1; $i++) {
		for ($j = $i+1; $j < $nifo; $j++) {
			$exclmatch=0;
			if (defined($excluded)){
				foreach (@excludedarray) {
					if (($_ eq $ARGV[$i].$ARGV[$j]) || $_ eq $ARGV[$j].$ARGV[$i]) {
						$exclmatch++;
					}

				}
			}

			unless ($exclmatch) {
                                $ifopair[$ij] = $ARGV[$i] . $ARGV[$j];
				$ij++;
			}
		}
	}

	


} else {
	die "Must specify interferometers";
}

$npairs=@ifopair;

if (defined($modifyfilter)){
        unless ($dorenormalize){
                die "Modify filter specified, but do renormalize not specified";        }else{
                @filter=split(/:/,$modifyfilter);
        }

}else{
        if ($dorenormalize){
                die "please specify a filter";
        }else{
                for($i=0;$i<$npairs;$i++){
                        push(@filter,0);
		}
        }
}


if (defined($excluded)){
#	print "npairs = $npairs, nexcl = $nexcl, nifo = $nifo\n";
	unless (($npairs +$nexcl) == $nifo * ($nifo-1) / 2 ){
		die "Must exclude available pairs";
	}
}

print "Creating dag file\n";

# set file names

$stochsub = "stochastic_multi.sub";
$postsub = "stochastic_multi_postproc.sub";
$combsub = "stochastic_multi_combine.sub";

# make directory for output of postproc
$postDir = "../output/${runid}/postProc/";
system("mkdir -p $postDir") && die "Can't create directory $postDir";

# make directory for err_ and out_ files
$outDir = "errout_${runid}/";
system("mkdir $outDir") && die "Can't create directory $outDir";


# get environment variables

$HOME = $ENV{"HOME"};
$LD_LIBRARY_PATH = $ENV{"LD_LIBRARY_PATH"};
$MATAPPS_TOP = $ENV{"MATAPPS_TOP"};

# set comment type
if ($shellcomments){
	$comment="#";
}else{
	$comment="%";
}

# open dag file

$filename = "stochastic_${runid}_condor.dag";

open( DAGFILE , ">$filename")
	or die "Could not open $filename";
	
print "=> $filename ...";


# dummy jobs to check params - one for each pair
for ( $i=0; $i < $npairs; $i++){
	print DAGFILE "JOB $i $stochsub\n";
	print DAGFILE "VARS $i paramsFile=\"../input/paramfiles/${ifopair[$i]}_${runid}_params.txt\" jobsFile=\"../input/jobfiles/${ifopair[$i]}_${runid}.txt\" jobNumber=\"$i\" home=\"${HOME}\" matapps_top=\"${MATAPPS_TOP}\" ld_library_path=\"${LD_LIBRARY_PATH}\" ifopair=\"$ifopair[$i]\" prevjobs=\"$i\" runid=\"$runid\"\n\n";
}

# first dummy job is parent, so that CTF archive extraction doesn't crash
# zeroth parent-child statement

if ($npairs>1){
	print DAGFILE "PARENT 0 CHILD ";
	for ($i=1;$i<$npairs;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "\n\n";
}
# loop over ifo pairs
$jobsprev=$npairs-1;

foreach $pair (@ifopair) {

# read job file and count njobs
	$jobfile = "../input/jobfiles/${pair}_${runid}.txt";

	unless (open (JOBFILE,"$jobfile")){
		die "Could not open $jobfile";
	}
	
	$njobs=0;

	while (<JOBFILE>){
		unless (substr($_,0,1) eq $comment){
			$njobs++;
		}
	}
# close job file
	close(JOBFILE);
	if ($njobs==0){
		die "no available jobs for pair $pair"
	}

#	print "$pair $njobs\n";
# for all jobs print an entry into dag file

	for ( $i=$jobsprev+1; $i <= $njobs+$jobsprev; $i++ ){
		print DAGFILE "JOB $i $stochsub\n";
		print DAGFILE "VARS $i paramsFile=\"../input/paramfiles/${pair}_${runid}_params.txt\" jobsFile=\"../input/jobfiles/${pair}_${runid}.txt\" jobNumber=\"${i}\" home=\"${HOME}\" matapps_top=\"${MATAPPS_TOP}\" ld_library_path=\"${LD_LIBRARY_PATH}\" ifopair=\"${pair}\" prevjobs=\"${jobsprev}\" runid=\"$runid\"\n\n";

	}
	$jobsprev+=$njobs;
}

# first parent - child relation
if ($npairs==1){
	print DAGFILE "PARENT 0 CHILD ";
	for ($i=1;$i<=$jobsprev;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "\n\n";
}else{
	print DAGFILE "PARENT ";
	for ($i=1;$i<$npairs;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "CHILD ";
	for ($i=$npairs;$i<=$jobsprev;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "\n\n";
}
# first postproc job done separately so that CTF extraction doesn't crash

# for all postproc jobs (loop over pairs) print an entry
$jobnum=$jobsprev;
$i=0;
foreach $pair (@ifopair){
	$jobnum++;
	print DAGFILE "JOB $jobnum $postsub\n";
	print DAGFILE "VARS $jobnum paramsFile=\"../input/paramfiles/${pair}_${runid}_params.txt\" jobsFile=\"../input/jobfiles/${pair}_${runid}.txt\" postDir=\"${postDir}\" dSigmaCut=\"$dsigmacut\" largeSigmaCutoff=\"$largesigmacut\" doRenormalize=\"$dorenormalize\" modifyFilter=\"\[$filter[$i]\]\" home=\"${HOME}\" matapps_top=\"${MATAPPS_TOP}\" ld_library_path=\"${LD_LIBRARY_PATH}\" ifopair=\"${pair}\" runid=\"$runid\"\n\n";
	$i++;
}

# second parent - child relation
if ($npairs==1){
	print DAGFILE "PARENT ";
        for ($i=1;$i<=$jobsprev;$i++){
                print DAGFILE "$i ";
        }
	$jobsp=$jobsprev+1;
	print DAGFILE "CHILD $jobsp\n\n";
}else{
	print DAGFILE "PARENT ";
	for ($i=2;$i<=$jobsprev;$i++){
		print DAGFILE "$i ";
	}
	$init=1+$jobsprev;
	print DAGFILE "CHILD $init\n\n";

	print DAGFILE "PARENT $init ";
	print DAGFILE "CHILD ";
	for ($i=2+$jobsprev;$i<=$jobnum;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "\n\n";
}

# final entry for combination job - only if > 1 pair
if ($npairs>1){
	$combjobnum=$jobnum+1;
	print DAGFILE "JOB $combjobnum $combsub\n";
	$ifostr=join("",@ARGV);
	print DAGFILE "VARS $combjobnum postDir=\"${postDir}\" ifos=\"${ifostr}\" home=\"${HOME}\" matapps_top=\"${MATAPPS_TOP}\" ld_library_path=\"${LD_LIBRARY_PATH}\"  runid=\"$runid\"\n\n";
	# third parent - child relation

	print DAGFILE "PARENT ";
	for ($i=1+$jobsprev;$i<=$jobnum;$i++){
		print DAGFILE "$i ";
	}
	print DAGFILE "CHILD $combjobnum\n\n";
}
# close file

close(DAGFILE);
print "done\n";
