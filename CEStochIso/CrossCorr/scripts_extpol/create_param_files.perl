#!/usr/bin/perl -w
# creates param files for multi ifo analysis 


use Getopt::Long;
#Getopt::Long::config('auto_abbrev');
Getopt::Long::config('no_ignore_case');

GetOptions('excluded|x:s' => \$excluded,				#x
	   'timeshifts|E:s' => \$timeshifts,				#E
	   'paramforce|p' => \$paramforce,				#p
	   'runid|i=s' => \$runid,					#i
	   'dofreqmask|m' => \$dofreqmask,				#m
	   'dohighpass|H' => \$dohighpass,				#H
	   'dooverlap|l' => \$dooverlap,				#l
	   'segmentduration|s:f' => \$segmentduration,			#s
	   'numsegmentsperinterval|n:i' => \$numsegmentsperinterval,	#n
	   'ignoremidsegment|I' => \$ignoremidsegment,			#I
	   'flow|f:f' => \$flow,					#f
	   'fhigh|F:f' => \$fhigh,					#F
	   'deltaF|D:f' => \$deltaF,					#D
	   'alphaexp|a:f' => \$alphaexp,				#a
	   'fref|R:f' => \$fref,					#R
	   'resamplerates|r:s' => \$resamplerates,			#r
	   'buffersecs|b:s' => \$buffersecs,				#b
	   'asqchannels|A:s' => \$asqchannels,				#A
	   'hannduration|S:s' => \$hannduration,			#S
	   'nresample|N:s' => \$nresample,				#N
	   'betaparam|B:s' => \$betaparam,				#B
	   'highpassfreq|P:s' => \$highpassfreq,			#P	
	   'highpassorder|O:s' => \$highpassorder,			#O
	   'freqstoremove|u:s' => \$freqstoremove,			#u
	   'nbinstoremove|U:s' => \$nbinstoremove,			#U
	   'shellcomments|C' => \$shellcomments,			#C
	   'frametypes|w:s' => \$frametypes,				#w
	   'framedurations|W:s' => \$framedurations,			#W
	   'injectionLocations|G:s' => \$injectionLocations,		#G
	   'injectionScales|J:s' => \$injectionScales,			#J
	   'signalType|q:s' => \$signalType,				#q
	   'simOmegaRef|Q:f' => \$simOmegaRef,				#Q
	   'doMonteCarlo|v' => \$doMonteCarlo,				#v
	   'doInjectionFromFile|V' => \$doInjectionFromFile,		#V
	   'numTrials|Z:i' => \$numTrials);				#Z

# Error catching and creation of pairs (excluding any if necessary)

if (defined($doMonteCarlo)) {
	if (defined($doInjectionFromFile)){
		die "Can't do injections from frame files and on-the-fly at the same time";
	}	
	unless (defined($simOmegaRef)) {
		die "doMontecarlo defined, but not simOmegaRef";
	}
	unless (defined($signalType)) {
		die "doMonteCarlo defined, but not signalType";
	}
        unless (defined($numTrials)) {
                die "doMonteCarlo defined, but not numTrials";
        }

} else {
	if (defined($simOmegaRef)) {
		die "simOmegaRef defined but not doMonteCarlo";
	}
	if (defined($signalType)) {
		die "signalType defined, but not doMondteCarlo";
	}
        if (defined($numTrials)) {
                die "numTrials defined, but not doMondteCarlo";
        }

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
					if ($_ eq $ARGV[$j].$ARGV[$i] || $_ eq $ARGV[$i].$ARGV[$j]) {
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

if (defined($timeshifts)){
        @tmp = split(/,/,$timeshifts);
        if (@tmp==$nifo){
                for($i=0;$i<$nifo;$i++){
                        $tshift{$ARGV[$i]} = $tmp[$i];
                }
        }elsif(@tmp==1){
                for($i=0;$i<$nifo;$i++){
                        $tshift{$ARGV[$i]} = $timeshifts;
                }
        }else{
                die "Must specify one timeshift per run or per ifo";
        }
}else{
        for($i=0;$i<$nifo;$i++){
                $tshift{$ARGV[$i]} = 0;
        }
}


if (defined($doInjectionFromFile)){
        if (defined($injectionLocations)){
                @tmp=split(/,/,$injectionLocations);
		if (@tmp==$nifo){
	                for($i=0;$i<$nifo;$i++){
               		        $injLoc{$ARGV[$i]} = $tmp[$i];
	                }
	        }elsif(@tmp==1){
			for($i=0;$i<$nifo;$i++){
	                        $injLoc{$ARGV[$i]} = $injectionLocations;
               		}
	        }else{
	                die "Number of injectionLocations must be unity or match number of ifos";
       		}
        }else{
                die "doInjectionFromFile defined, but not injectionLocations";
        }
        if (defined($injectionScales)){
                @tmp=split(/,/,$injectionScales);
		if (@tmp==$nifo){                         
			for($i=0;$i<$nifo;$i++){
                                $injScale{$ARGV[$i]} = $tmp[$i];    
			}
                }elsif(@tmp==1){
                        for($i=0;$i<$nifo;$i++){
                                $injScale{$ARGV[$i]} = $injectionScales;
                        }
                }else{
                        die "Number of injectionLocations must be unity or match number of ifos";
                }
	}else{
                die "doInjectionFromFile defined, but not injectionScales";
        }

}else{
        if (defined($injectionScales)){
                die "injectionScales defined, but not doInjectionFromFile";
        }
        if (defined($injectionLocations)){
                die "injectionLocations defined, but not doInjectionFromFile";
        }
}


unless(defined($alphaexp)){
	$alphaexp=0;
}

unless(defined($fref)){
	$fref=100;
}

unless(defined($flow)){
	$flow=40;
}

unless(defined($fhigh)){
	$fhigh=500;
}

unless(defined($deltaF)){
	$deltaF=0.25;
}

unless (defined($numsegmentsperinterval)){
	$numsegmentsperinterval=1;
}

unless (defined($segmentduration)){
	$segmentduration=60;
}

#if (defined($dohighpass1)){
#	unless (defined($highpassfreq1)){
#		$highpassfreq1=40;
#	}
#	unless (defined($highpassorder1)){
#		$highpassorder1=6;
#	}
#}

#if (defined($dohighpass2)){
#	unless (defined($highpassfreq2)){
#		$highpassfreq2=40;
#	}
#	unless (defined($highpassorder2)){
#		$highpassorder2=6;
#	}
#}

unless (defined($runid)) {
	die "Please give a unique identifier for your analysis";
}


if (defined($excluded)){
	#print "npairs = $npairs, nexcl = $nexcl, nifo = $nifo\n";
	#print "@ifopair/n";	
	unless (($npairs +$nexcl) == $nifo * ($nifo-1) / 2 ){
		die "Must exclude available pairs";
	}
}

if (defined($resamplerates)){
	@tmp = split(/,/,$resamplerates);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$resample{$ARGV[$i]} = $tmp[$i];
		}
	}else{
		die "Number of resample rates must match number of ifos";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		if($ARGV[$i] eq 'V1'){
			$resample{$ARGV[$i]} = 1000;
		}else{
			$resample{$ARGV[$i]} = 1024;
		}
	}
}
	
if (defined($buffersecs)){
	@tmp = split(/,/,$buffersecs);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$buffer{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$buffer{$ARGV[$i]} = $buffersecs;
		}
	}else{
		die "Number of buffers must be unity or match number of ifos";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$buffer{$ARGV[$i]} = 1;
	}
}

if (defined($asqchannels)){
	@tmp = split(/,/,$asqchannels);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$asq{$ARGV[$i]} = $tmp[$i];
		}
        }elsif(@tmp==1){
                for($i=0;$i<$nifo;$i++){
                        $asq{$ARGV[$i]} = $asqchannels;
                }
	}else{
		die "Number of asq channels must match number of ifos or be unity";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		if($ARGV[$i] eq 'V1'){
			$asq{$ARGV[$i]} = "h_20000Hz";
		}elsif($ARGV[$i] eq 'G1'){
			$asq{$ARGV[$i]} = "DER_DATA_H";
		}else{
			$asq{$ARGV[$i]} = "LSC-STRAIN";
		}
	}
}

if (defined($frametypes)){
	@tmp = split(/,/,$frametypes);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$frtype{$ARGV[$i]} = $tmp[$i];
		}
        }elsif(@tmp==1){
                for($i=0;$i<$nifo;$i++){
                        $frtype{$ARGV[$i]} = $frametypes;
                }
	}else{
		die "Number of frame types must match number of ifos";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		if($ARGV[$i] eq 'V1'){
			$frtype{$ARGV[$i]} = "HrecV2_20000Hz";
		}elsif($ARGV[$i] eq 'G1'){
			$frtype{$ARGV[$i]} = "RDS_C01_L3";
		}else{
			$frtype{$ARGV[$i]} = "$ARGV[$i]_RDS_C03_L2";
		}
	}
}
if (defined($framedurations)){
	@tmp = split(/,/,$framedurations);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$frdurn{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$frdurn{$ARGV[$i]} = $framedurations;
		}
	}else{
		die "Must specify one frame duration per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$frdurn{$ARGV[$i]} = -1;
	}
}


if (defined($hannduration)){
	@tmp = split(/,/,$hannduration);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$hann{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$hann{$ARGV[$i]} = $hannduration;
		}
	}else{
		die "Must specify one hann duration per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$hann{$ARGV[$i]} = $segmentduration;
	}
}

if (defined($nresample)){
	@tmp = split(/,/,$nresample);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$nres{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$nres{$ARGV[$i]} = $nresample;
		}
	}else{
		die "Must specify one nresample per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$nres{$ARGV[$i]} = 10;
	}
}

if (defined($betaparam)){
	@tmp = split(/,/,$betaparam);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$beta{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$beta{$ARGV[$i]} = $betaparam;
		}
	}else{
		die "Must specify one betaparam per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$beta{$ARGV[$i]} = 5;
	}
}

if (defined($highpassfreq)){
	@tmp = split(/,/,$highpassfreq);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$hpfreq{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$hpfreq{$ARGV[$i]} = $highpassfreq;
		}
	}else{
		die "Must specify one highpassfreq per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$hpfreq{$ARGV[$i]} = 40;
	}
}

if (defined($highpassorder)){
	@tmp = split(/,/,$highpassorder);
	if (@tmp==$nifo){
		for($i=0;$i<$nifo;$i++){
			$hporder{$ARGV[$i]} = $tmp[$i];
		}
	}elsif(@tmp==1){
		for($i=0;$i<$nifo;$i++){
			$hporder{$ARGV[$i]} = $highpassorder;
		}
	}else{
		die "Must specify one highpassorder per run or per ifo";
	}
}else{
	for($i=0;$i<$nifo;$i++){
		$hporder{$ARGV[$i]} = 6;
	}
}

if(defined($dofreqmask)){
	if(defined($nbinstoremove)){
		unless(defined($freqstoremove)){
			die "Must define frequencies to remove if bins to remove are defined";
		}
	}

	if (defined($freqstoremove)){
		@tmp = split(/:/,$freqstoremove);
		if (@tmp==$nifo){
			for($i=0;$i<$nifo;$i++){
				$freqrem{$ARGV[$i]} = $tmp[$i];
			}
		}elsif(@tmp==1){
			for($i=0;$i<$nifo;$i++){
				$freqrem{$ARGV[$i]} = $freqstoremove;
			}
		}else{
			die "Must specify one freqstoremove list per run or per ifo";
		}
	}else{
		die "Must specify list of frequencies to remove";
	}	

	if (defined($nbinstoremove)){
		@tmp = split(/:/,$nbinstoremove);
		if (@tmp==$nifo){
			for($i=0;$i<$nifo;$i++){
				$nbinrem{$ARGV[$i]} = $tmp[$i];
				@t1=split(/,/ , $nbinrem{$ARGV[$i]});
				$l1=@t1;
				@t2=split(/,/ , $freqrem{$ARGV[$i]});
				$l2=@t2;
				unless ($l1==$l2){
					die "Length of nbins to remove must match length of freqs to remove";
				}
			}
		}elsif(@tmp==1){
			for($i=0;$i<$nifo;$i++){
				$nbinrem{$ARGV[$i]} = $nbinstoremove;
				@t1=split(/,/ , $nbinrem{$ARGV[$i]});
				$l1=@t1;
				@t2=split(/,/ , $freqrem{$ARGV[$i]});
				$l2=@t2;
				unless ($l1==$l2){
					die "Length of nbins to remove must match length of freqs to remove";
				}
			
			}
		}else{
			die "Must specify one nbinstoremove list per run or per ifo";
		}

	}else{
		if(defined($freqstoremove)){
			for($i=0;$i<$nifo;$i++){
				$_=$freqrem{$ARGV[$i]};
				s/\d{1,}/1/g;
				$nbinrem{$ARGV[$i]}=$_;
			
			}
		}
	}
}

# System variables
$user=`echo \$USER`;
chomp($user);
$hostname=`echo \$HOSTNAME`;
chomp($hostname);
$now=`date`;
chomp($now);

# If required, create paramfiles. Default is to create file if one doesn't exist
# and not to if it does. This can be forced to overwrite existing files.

$paramdir="../input/paramfiles/";
print "creating param files:\n";	

for ($i=0;$i<$npairs;$i++){

	$ifo1=substr($ifopair[$i],0,2);
	$ifo2=substr($ifopair[$i],2,2);

	$paramfilename= "$paramdir$ifopair[$i]_${runid}_params.txt";
	print "=> $paramfilename ...";		
	if ($paramforce || !(-e $paramfilename)) {
		open( PARAMFILE , ">$paramfilename")
			or die "Could not open $paramfilename";

		print PARAMFILE "% param file generated by running\n";
		print PARAMFILE "% $0 @ARGV\n";
		print PARAMFILE "% Created by $user on $hostname at $now\n"; 
		print PARAMFILE "\n% flags for optional operations\n";
		if ($dofreqmask){
			print PARAMFILE "doFreqMask true\n";
		}else{
			print PARAMFILE "doFreqMask false\n";
		}
		if ($dohighpass){
			print PARAMFILE "doHighPass1 true\n";
#			print PARAMFILE "highPassFreq1 $highpassfreq1\n";
#			print PARAMFILE "highPassOrder1 $highpassorder1\n";
		}else{
			print PARAMFILE "doHighPass1 false\n";
		}
		if ($dohighpass){
			print PARAMFILE "doHighPass2 true\n";
#			print PARAMFILE "highPassFreq2 $highpassfreq2\n";
#			print PARAMFILE "highPassOrder2 $highpassorder2\n";
		}else{
			print PARAMFILE "doHighPass2 false\n";
		}
		if ($doMonteCarlo){
			print PARAMFILE "doMonteCarlo true\n"; 
		}else{
                        print PARAMFILE "doMonteCarlo false\n";
		}
		if ($doInjectionFromFile){
                        print PARAMFILE "doInjFromFile1 true\n";
			print PARAMFILE "doInjFromFile2 true\n";
		}else{
                        print PARAMFILE "doInjFromFile1 false\n";
                        print PARAMFILE "doInjFromFile2 false\n";
		}	
		if ($dooverlap){
			print PARAMFILE "doOverlap true\n";
		}else{
			print PARAMFILE "doOverlap false\n";
		}
		print PARAMFILE "heterodyned false\n";
		print PARAMFILE "suppressFrWarnings true\n";
		print PARAMFILE "writeNaiveSigmasToFiles true\n";
		print PARAMFILE "writeResultsToScreen true\n";
		print PARAMFILE "writeStatsToFiles true\n";
		print PARAMFILE "writeSpectraToFiles true\n";
		print PARAMFILE "writeSensIntsToFiles true\n";
		print PARAMFILE "writeOptimalFiltersToFiles true\n";
		print PARAMFILE "writeCalPSD1sToFiles true\n";
		print PARAMFILE "writeCalPSD2sToFiles true\n";
		print PARAMFILE "\n% ifo names\n";	
		print PARAMFILE "ifo1 $ifo1\n";
		print PARAMFILE "ifo2 $ifo2\n";
		print PARAMFILE "\n% segment duration (sec)\n";
		print PARAMFILE "segmentDuration $segmentduration\n";
		print PARAMFILE "\n%timeshifts\n";
		if ($tshift{$ifo1}==0){
			print PARAMFILE "doShift1 false\n";
		}else{
			print PARAMFILE "doShift1 true\n";
			print PARAMFILE "ShiftTime1 $tshift{$ifo1}\n";
		}
                if ($tshift{$ifo2}==0){
                        print PARAMFILE "doShift2 false\n";
                }else{
                        print PARAMFILE "doShift2 true\n";
                        print PARAMFILE "ShiftTime2 $tshift{$ifo2}\n";
                }

		print PARAMFILE "\n% parameters for sliding psd estimation:\n";
		print PARAMFILE "% numSegmentsPerInterval should be odd; ignoreMidSegment is a flag\n";
		print PARAMFILE "% that allows you to ignore (if true) or include (if false) the\n";
		print PARAMFILE "% analysis segment when estimating power spectra\n";
		print PARAMFILE "numSegmentsPerInterval $numsegmentsperinterval\n";
		if ($ignoremidsegment){
			print PARAMFILE "ignoreMidSegment true\n";
		}else{
			print PARAMFILE "ignoreMidSegment false\n";
		}		
		print PARAMFILE "\n% freq resolution and freq cutoffs for CC statistic sum (Hz)\n";
		print PARAMFILE "flow $flow\n";
		print PARAMFILE "fhigh $fhigh\n";
		print PARAMFILE "deltaF $deltaF\n";
		print PARAMFILE "\n% params for Omega_gw (power-law exponent and reference freq in Hz)\n";
		print PARAMFILE "alphaExp $alphaexp\n";
		print PARAMFILE "fRef $fref\n";
		print PARAMFILE "\n% resample rate (Hz)\n";
		print PARAMFILE "resampleRate1 $resample{$ifo1}\n";
		print PARAMFILE "resampleRate2 $resample{$ifo2}\n";		
		print PARAMFILE "\n% buffer added to beginning and end of data segment to account for\n";
		print PARAMFILE "% filter transients (sec)\n";
		print PARAMFILE "bufferSecs1 $buffer{$ifo1}\n";
		print PARAMFILE "bufferSecs2 $buffer{$ifo2}\n";
		print PARAMFILE "\n% ASQ channel\n";
		print PARAMFILE "ASQchannel1 $asq{$ifo1}\n";
		print PARAMFILE "ASQchannel2 $asq{$ifo2}\n";
		if (defined($doInjectionFromFile)){
			print PARAMFILE "injChannel1 strain\n";
                        print PARAMFILE "injChannel2 strain\n";
		}
		print PARAMFILE "\n% frame type and duration\n";
		print PARAMFILE "frameType1 $frtype{$ifo1}\n";
		print PARAMFILE "frameType2 $frtype{$ifo2}\n";
		print PARAMFILE "frameDuration1 $frdurn{$ifo1}\n";
		print PARAMFILE "frameDuration2 $frdurn{$ifo2}\n";
		print PARAMFILE "\n% duration of hann portion of tukey window \n";
		print PARAMFILE "\n% (hannDuration = segmentDuration is a pure hann window)\n";
		print PARAMFILE "hannDuration1 $hann{$ifo1}\n";
		print PARAMFILE "hannDuration2 $hann{$ifo2}\n";
		print PARAMFILE "\n% params for matlab resample routine\n";
		print PARAMFILE "nResample1 $nres{$ifo1}\n";
		print PARAMFILE "nResample2 $nres{$ifo2}\n";
		print PARAMFILE "betaParam1 $beta{$ifo1}\n";
		print PARAMFILE "betaParam2 $beta{$ifo2}\n";
		print PARAMFILE "\n% params for high-pass filtering (3db freq in Hz, and filter order) \n";
		print PARAMFILE "highPassFreq1 $hpfreq{$ifo1}\n";
		print PARAMFILE "highPassFreq2 $hpfreq{$ifo2}\n";
		print PARAMFILE "highPassOrder1 $hporder{$ifo1}\n";
		print PARAMFILE "highPassOrder2 $hporder{$ifo2}\n";
		print PARAMFILE "\n% coherent freqs and number of freq bins to remove if doFreqMask=true;\n";
		print PARAMFILE "% NOTE: if an nBin=0, then no bins are removed even if doFreqMask=true\n";
		print PARAMFILE "% (coherent freqs are typically harmonics of the power line freq 60Hz\n";
		print PARAMFILE "% and the DAQ rate 16Hz)\n";
		if(defined($dofreqmask)){
			@frem1=split(/,/,$freqrem{$ifo1});
			@frem2=split(/,/,$freqrem{$ifo2});
			@nrem1=split(/,/,$nbinrem{$ifo1});	
			@nrem2=split(/,/,$nbinrem{$ifo2});	
			$k=0;
			$j=0;
			foreach $fr1 (@frem1){
				foreach $fr2 (@frem2){
					if ($fr1 == $fr2){
					#	$indx[$j]=$k;
					#	$j++;
						if ($nrem1[$k]>$nrem2[$k]){
                                                        $frem2[$k]=-1;
                                                        $nrem2[$k]=-1;
						}else{
							$frem1[$k]=-1;
							$nrem1[$k]=-1;
						}
					}
				}
				$k++;
			
			}
		
			$fr=stripneg1(@frem1);
			$nr=stripneg1(@nrem1);
			$fr2=join(",",$fr,$freqrem{$ifo2});
			$nr2=join(",",$nr,$nbinrem{$ifo2});
			print PARAMFILE "freqsToRemove $fr2\n";
			print PARAMFILE "nBinsToRemove $nr2\n";
		}
		print PARAMFILE "\n% number of trials for monte carlo simulations (if doMonteCarlo = true)\n";
		if ($doMonteCarlo) {
			print PARAMFILE "numTrials $numTrials\n";
		}else{
                        print PARAMFILE "numTrials -1\n";
		}
		print PARAMFILE "\n% type of SB signal to simulate \n";
		print PARAMFILE "% (const for Omega_gw=const, white for Omega_gw propto f^3)\n";
		if ($doMonteCarlo){
			print PARAMFILE "signalType $signalType\n";
		}else{
                        print PARAMFILE "signalType const\n";
		}
		
		print PARAMFILE "\n% value of Omega_gw(f_Ref) for simulated SB signal\n";
		if ($doMonteCarlo){
			print PARAMFILE "simOmegaRef $simOmegaRef\n";
		}else{
                        print PARAMFILE "simOmegaRef 1\n";
		}
		
		if ($doInjectionFromFile){
			print PARAMFILE "\n% strain scale factor for injected signals\n";
			print PARAMFILE "injScale1 $injScale{$ifo1}\n"; 
                        print PARAMFILE "injScale2 $injScale{$ifo2}\n";
                        #print "injScale1 $injScale{$ifo1}\n";
                        #print "injScale2 $injScale{$ifo2}\n";

		}	
		print PARAMFILE "\n% prefix convention for comments in job files\n";
		if ($shellcomments){
			print PARAMFILE "jobsFileCommentStyle shell\n";
		}else{
			print PARAMFILE "jobsFileCommentStyle matlab\n";
		}

		print PARAMFILE "\n% calibration filenames\n";
		print PARAMFILE "alphaBetaFile1 none\n";
		print PARAMFILE "alphaBetaFile2 none\n"; 
		print PARAMFILE "calCavGainFile1 none\n"; 
		print PARAMFILE "calCavGainFile2 none\n";
		print PARAMFILE "calResponseFile1 none\n";
		print PARAMFILE "calResponseFile2 none\n";
	
		print PARAMFILE "\n% path to cache files\n";
		system("mkdir -p ../input/cachefiles/$runid/$ifopair[$i]") && die "couldn't create directory ../input/cachefiles/$runid/$ifopair[$i]";

		print PARAMFILE "gpsTimesPath1 ../input/cachefiles/$runid/$ifopair[$i]/\n";
		print PARAMFILE "gpsTimesPath2 ../input/cachefiles/$runid/$ifopair[$i]/\n";
		print PARAMFILE "frameCachePath1 ../input/cachefiles/$runid/$ifopair[$i]/\n";
		print PARAMFILE "frameCachePath2 ../input/cachefiles/$runid/$ifopair[$i]/\n";
	
		if ($doInjectionFromFile){
			print PARAMFILE "\n% path to cache files for injection\n";
			print PARAMFILE "injGPSTimesPath1 <auto>\n";
			print PARAMFILE "injGPSTimesPath2 <auto>\n";
			system("mkdir -p ../input/cachefiles/${runid}_inj/$ifopair[$i]/");
			if (defined($shellcomments)){
				system("./create_inj_cache.perl -i ${runid} -p ${ifopair[$i]} -L $injLoc{$ifo1} -C");
			}else{
				system("./create_inj_cache.perl -i ${runid} -p ${ifopair[$i]} -L $injLoc{$ifo2}");
			}
			print PARAMFILE "injFrameCachePath1 ../input/cachefiles/${runid}_inj/$ifopair[$i]/\n";
			print PARAMFILE "injFrameCachePath2 ../input/cachefiles/${runid}_inj/$ifopair[$i]/\n";

		}
	
		print PARAMFILE "\n% prefix for output filename\n";
		system("mkdir -p ../output/$runid/$ifopair[$i]/") && die "Can't create directory ../output/$runid/$ifopair[$i]/";
		print PARAMFILE "outputFilePrefix ../output/$runid/$ifopair[$i]/$ifopair[$i]_$runid\n";
		print "done\n";

	}else{
		print "file exists - skipping\n";
	}	

}


sub stripneg1 {
	$_=join(",",@_);
	s/,-1,/,/g;
	s/,-1//g;
	s/-1,//g;
	if ($_ eq "-1"){
		$_="";
	}
	return $_;
}
