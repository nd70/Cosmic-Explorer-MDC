#! /usr/bin/perl
# E. Thrane
# combines jobfiles in jobfiles/ to create jobfile_master
#

$i=1;
# Find the total number of jobfiles
$jmax=`ls ./jobfiles/ | sed 's/jobfile_//g' | sort -bn | tail -1`;
open(master,">jobfile_master");
for ($j=1; $j<=$jmax; $j++) {
    $jobfile = "./jobfiles/jobfile_$j";
    open(jobfile,$jobfile);
    while (<jobfile>) {
	chomp($_);
	$start=$_;  $start=~s/.* (.*) (.*) (.*)/\1/;
	$stop=$_;  $stop=~s/.* (.*) (.*) (.*)/\2/;
	$dur=$_;  $dur=~s/.* (.*) (.*) (.*)/\3/;
	print master "$i $start $stop $dur\n";
	$i = $i + 1;
    }
    close(jobfile)
}
close(master);
