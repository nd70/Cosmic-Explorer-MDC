#!/usr/bin/perl -w
#
# $Id: buildJobFile_mat.perl,v 1.1 2009-09-08 12:49:47 elr Exp $
#

use Fcntl;
use POSIX qw(log10 ceil);
use Carp;
require 'dumpvar.pl';
use Getopt::Long;
Getopt::Long::config('auto_abbrev');

GetOptions('start=i' => \$gpsstarttime,
	   'end=i' => \$gpsstoptime,
	   'min=i'=> \$minlength,
	   'pgplus:i'=> \$pgplus,
	   'npgplus:i'=> \$npgplus,
	   'max=i'=> \$maxlength);

$pglength = 600;
$pgperiod = 6370;
$pgstart = 729273613  % $pgperiod;

if (defined $pgplus) {
    if (defined $npgplus) {
	die "Can't specify both playground and non-playground";
    } else {
	$pgsellength = $pglength + 2*$pgplus;
	$pgselstart = ($pgstart-$pgplus) % $pgperiod;
    }
} elsif (defined $npgplus) {    
    $pgsellength = $pgperiod - $pglength + 2*$npgplus;
    $pgselstart = ($pgstart+$pglength-$npgplus) % $pgperiod;
}

if (defined $pgsellength) {
    # warn $pgsellength;
    # warn $pgselstart;
    # warn $minlength;
    if ($pgsellength <= 0) {
	die "playground selection has nonpositive length $pgsellength";
    } elsif (defined $minlength and $pgsellength < $minlength) {
	die "playground selection length $pgsellength less than min length $minlength";
    }
}

if (defined $maxlength and defined $minlength and $maxlength <= $minlength){
	die "maximum job length $maxlength must be greater than minimum length $minlength";
}

print "% Segment list produced by running\n";
print "% $0 @ARGV\n";
print "% On segwizard output which was\n";
while(defined ($_ = <STDIN>)) {
    if (substr($_,0,1) eq '#') {
	if (substr($_,0,5) eq '# Don') {
	    if (defined $gpsstarttime or defined $gpsstoptime) {
		print "% Restrict to GPS time range $gpsstarttime - $gpsstoptime\n";
	    }
	    if (defined $pgplus) {
		if ($pgplus == 0) {
		    print "% Use only playground data\n";
		} else {
		    print "% Use only playground data, plus $pgplus sec on either side\n";
		}
	    }
	    elsif (defined $npgplus) {
		if ($npgplus == 0) {
		    print "% Use only non-playground data\n";
		} else {
		    print "% Use only non-playground data, plus $npgplus sec on either side\n";
		}
	    }
	    if (defined $minlength) {
		print "% Minimum segment length to analyze: $minlength\n";
	    }
	    if (defined $maxlength) {
		print "% Maximum segment length to analyze: $minlength\n";
            }
	    $segstring = '';
	    $totjoblength = 0;
	}
	if (substr($_,0,5) eq '# Tot') {
	     # have to wait for last two lines of header until know total time
	} else {
	    	print "%";
		print $_;
	}
    } else {
	($segnum,$jobstart,$jobstop,$joblength) = split;
	if ($joblength != $jobstop - $jobstart) {
	    die "Segment length $joblength disagrees with stop $jobstop minus start $jobstart";
	}
	if (defined $gpsstarttime) {
	    next if $gpsstarttime >= $jobstop;
	    if ($gpsstarttime > $jobstart) {
		$jobstart = $gpsstarttime;
		$joblength = $jobstop - $jobstart;
	    }
	    undef $gpsstarttime; # don't need any more
	}
	if (defined $gpsstoptime) {
	    last if ($gpsstoptime <= $jobstart);
	    if ($gpsstoptime < $jobstop) {
		$jobstop = $gpsstoptime;
		$joblength = $jobstop - $jobstart;		
	    }
	}
	# warn $joblength;
	next if $joblength < $minlength;
	if (defined $maxlength) {
		$subjobstart = $jobstart;
		$subjobstop = $jobstart+$maxlength;

		while ($subjobstop <= $jobstop){
			if (defined $pgsellength) {
	                    &step_through_pg_sels($subjobstart,$subjobstop);
			} else {
	                    $segstring .=
        	                sprintf("%4d  %9d  %9d %6d\n",
                	                $segnum,$subjobstart,$subjobstop,$maxlength);
	                    $totjoblength += $maxlength;
        	        }
			$subjobstart += $maxlength;
			$subjobstop += $maxlength;
		}
		
		if ($subjobstop > $jobstop) {
			$subjobstop = $jobstop;
			$subjoblength = $subjobstop - $subjobstart;
			if (defined $pgsellength) {
                            &step_through_pg_sels($subjobstart,$subjobstop);
                        } else {                             $segstring .=
                                sprintf("%4d  %9d  %9d %6d\n",
                                        $segnum,$subjobstart,$subjobstop,$subjoblength);
                            $totjoblength += $subjoblength;
                        }

		}

	} else {
		if (defined $pgsellength) {	    
		    &step_through_pg_sels($jobstart,$jobstop);
		} elsif (not defined $minlength or  $joblength >= $minlength) {
		    $segstring .= 
			sprintf("%4d  %9d  %9d %6d\n",
				$segnum,$jobstart,$jobstop,$joblength);
		    $totjoblength += $joblength;
		}
	}
    }
}
print "% Total duration of segments to be analyzed: $totjoblength\n";
print "%seg    start      stop    duration\n";
print $segstring;

sub step_through_pg_sels {
    $ti = shift @_;
    $tf = shift @_;

    $tjump = ($pgselstart-$ti) % $pgperiod;
    # warn $tjump;
    $thislength = $pgsellength - ($pgperiod - $tjump);
    if ($ti+$thislength > $tf) {
	$thislength = $tf-$ti;
    }
    # warn $thislength;
    if ($thislength >= $minlength) {
	$segstring .= 
	    sprintf("%4d  %9d  %9d %6d\n",
		    $segnum,$ti,$ti+$thislength,$thislength);	
	$totjoblength += $thislength;
    }
    $ti += $tjump;
    while ($ti < $tf) {
	$thisend = $ti + $pgsellength;
	# warn "$ti $thisend $tf";
	if ($thisend <= $tf) {
	    $segstring .= 
		sprintf("%4d  %9d  %9d %6d\n",
			$segnum,$ti,$thisend,$pgsellength);
	    $totjoblength += $pgsellength;
	} else {
	    $thislength = $tf - $ti;
	    if ($thislength >= $minlength) {
		$segstring .= 
		    sprintf("%4d  %9d  %9d %6d\n",
			    $segnum,$ti,$tf,$thislength);
		$totjoblength += $thislength;
	    }
	}
	$ti += $pgperiod;
    }
}
