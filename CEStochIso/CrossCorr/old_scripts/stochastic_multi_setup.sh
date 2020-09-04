#!/bin/bash
# runs all the scripts required to set up and submit the stochastic multi ifo
# analysis. Matapps must be installed, and compiled in the stochastic subdir.

# Get options

#set defaults
PARAM="1"
CACHE="1"
JOB="1"
DAG="1"
SUBMIT="1"
SHELL="0"
PG="0"
FORCEDAG="0"
FORCEJOB="0"
FORCECACHE="0"
FORCEPARAM="0"
RENORM="0"
PARAMOPT=""
CACHEOPT=""
JOBOPT=""
DAGOPT=""

usage="help"

while getopts  ":a:b:cdef:ghi:jklmn:opq:r:s:t:u:vw:x:yzA:B:CD:E:F:G:HIJ:K:L:M:N:O:P:Q:R:S:T:U:VW:XYZ:" flag
do
#  echo "$flag" $OPTIND $OPTARG
  case "$flag" in
  	b) buffer=$OPTARG
           PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
	c) FORCECACHE="1"
	   CACHEOPT="-f ${CACHEOPT}";;
	d) FORCEDAG="1";;
	e) CACHE="0";;
	g) DAG="0";;
	h) echo "Usage stochastic_multi_setup [options], where options are: (optional are in square brackets)

[a]	exponent of power law, default=0
[b]	buffer seconds, one value, or comma separated list (one per ifo). If not specified, default=1
[c]	force creation of cache files, even if they exist
[d]	force creation of dag files, even if they exist 
[e]	don't create cache files
f	minimum frequency to be analysed
[g]	don't create dag
h	show help
i	unique identifier for this analysis
[j]	force creation of job files, even if they exist
[k]	don't create job files
[l]	overlap time segments
[m]	mask frequencies specified by -u/-U
[n]	number of segments per interval to be used to estimate PSDs (if not specified, default=1)
[o]	don't create param files
[p]	force creation of param files, even if they exist
[q]	type of signal to simulate, white or flat
[r]	resample rates for each detector, one value or comma separated list. If not specified default is 1000 for Virgo and 1024 for LIGO and GEO
[s]	duration of analysis segments (if not specified, default=60)
t	start of analysis
[u]	frequency bins to remove when doing frequency mask. One comma separated list, or one list per IFO, separated by colons, eg. 1,1,1:2,3,4:1,2,1
[v]	do monte carlo simulation for SW injection
[w]	frameTypes to be analysed, one value or comma separated list. If not specified, default is HrecV2_20000Hz for Virgo, RDS_C01_L3 for GEO and xx_RDS_C03_L2 for H1, L1 or H2
[x]	pairs not to be analysed, comma separated list, eg H1H2,G1L1
[y]	only choose playground segments
[z]	renormalise the data, using the modify filters in files specified by the prefix -M
[A]	ASQ channels to be analysed, one value or comma separated list. If not specified, default is h_20000Hz for Virgo, DER_DATA_H for GEO and LSC-STRAIN for LIGO
[B]	beta parameter for resampleing routine, single value or comma separated list. If not specified, default is 5
[C]	comment style in jobfile is shell (i.e. #), default is 'matlab' (i.e. %). This is useful when one is using job files written by other scripts
[D]	resolution of optimal filter. If not defined, default is 0.25
[E]	timeshift for each detector (comma separated list)
[F]	maximum frequency
[G]	location of frame files containing injections
[H]     high-pass filter the data (default is to do highpass)
[I]	ignore the centre segment of the interveal when finding PSDs
[J]	scale for SW injectino from file, one value or comma separated list
[K]	the delta sigma cut to be applied in postprocessing
[L]	the large sigma cut to be applied in postprocessing
[M]	the sufffix of the files in which the modify filter is stored
[N]	the nresample parameter for the resampling routine. If undefined, default=10
[O]	the order of the high pass filtering, single value or comma separated list. If undefined, default=6
[P]	the 3db frequency for high pass filtering, single value or comma separated list. If undefined, default=40
[Q]	Omega_R for on-the-fly SW injections
[R]	reference frequnecy. Default=100
[S]	duration of hann portion of window s=S is pure hann window. Default is to use value given for segment duration
[T]	end of the analysis
[U]	the number of bins to remove for each list of frequencies to remove
[V]	do a SW injection using frame files specified by -G
[W]	durations of frame files, signal value or comma-separeted list. Can specify "-1" if the frame durations are to be found from the frame file names. Default=-1
[X]	don't submit the dags, just set up the analysis
[Y]	only do non-playground segments
[Z]	number of trials for the on-the-fly injections";;
	i) runid=$OPTARG
	   DAGOPT="-$flag $OPTARG ${DAGOPT}"
           PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
	k) JOB="0";;
	j) FORCEJOB="1";;
	n) nspi=$OPTARG
	   PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
        o) PARAM="0";;
	p) FORCEPARAM="1"
	   PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
	s) duration=$OPTARG
           PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
	t) JOBOPT="--start $OPTARG $JOBOPT"
	   gpsStartTime=$OPTARG;;
	v) PARAMOPT="-$flag $OPTARG ${PARAMOPT}"
	   JOBOPT="--max 242 $JOBOPT";;
	x) excluded=$OPTARG
	   DAGOPT="-$flag $OPTARG ${DAGOPT}"
           PARAMOPT="-$flag $OPTARG ${PARAMOPT}";;
	y) PG="1";;
	z) RENORM="1"
	   DAGOPT="-$flag $OPTARG ${DAGOPT}";;
	C) SHELL="1"
	   PARAMOPT="-$flag $OPTARG ${PARAMOPT}"
	   DAGOPT="-$flag $OPTARG ${DAGOPT}";;
	K) DAGOPT="-$flag $OPTARG ${DAGOPT}";;
	L) DAGOPT="-$flag $OPTARG ${DAGOPT}";;
	M) modifyFilterSuffix=$OPTARG;;
	T) JOBOPT="--end $OPTARG $JOBOPT"
	   gpsEndTime=$OPTARG;;
	X) SUBMIT="0";;
	Y) PG="-1";;
	\?) echo "Invalid option: -$OPTARG"
#"$usage 1 $flag $OPTARG"
	    exit 1;;
  	*) if [ ! `echo afhilmnpqrsuwxABDEFGHIJNOPQRSUVWZ | grep -q $flag` ]; then 
	     PARAMOPT="-$flag $OPTARG ${PARAMOPT}"
	   else
	     echo "$usage 2 $flag $OPTARG"
	     exit 1
	   fi
	   ;;

  esac
done
# Check $MATAPPS_TOP is set
if [ -z "$MATAPPS_TOP" ]; then
	echo "Please set the \$MATAPPS_TOP environment variable"
	exit 1
fi

# Check that the analysis has a uinique identifier
if [ -z "$runid" ]; then
	echo "Please enter a unique identifier for your analysis"
	exit 1
fi

# shift args so that $@/$* contain just the ifos
snum=$(( $OPTIND - 1 ))
shift $snum

# Make sure we have some ifos to use
if [ "$#" -lt 1 ]; then
	echo "Please provide the interferometers you wish to use"
	exit 1
fi

# Check for conflicts between forcing files and no files
if [ $FORCEDAG -gt $DAG ]; then
	echo "Conflict between forcedag and nodag"
	exit 1
fi
if [ $FORCECACHE -gt $CACHE ]; then
	echo "Conflict between forcecache and nocache"
	exit 1
fi
if [ $FORCEJOB -gt $JOB ]; then
	echo "Conflict between forcejob and nojob"
	exit 1
fi
if [ $FORCEJOB -gt $JOB ]; then
	echo "Conflict between forcejob and nojob"
	exit 1
fi
if [ $FORCEPARAM -gt $PARAM ]; then
	echo "Conflict between forceparam and noparam"
	exit 1
fi

# check that modifyfilter file is given, otherwise set to default
if [ -z "$modifyFilterSuffix" ]; then
	if [ "$RENORM" == "1" ]; then
		modifyFilterSuffix="_${runid}_filter.txt"
	fi

else
	if [ "$RENORM" == "0" ]; then
		echo Modify filter file given, but not set to renormalize
		exit 1
	fi
fi

# set minimum value for length of data (using defaults where necessary)
# (numSegsPerInterval*segment length + 2*buffer)
if [ -z "$buffer" ]; then
	buffer="1"
fi
if [ -z "$duration" ]; then
	duration="60"
fi
if [ -z "$nspi" ]; then
	nspi="3"
fi
nspi=$(( $nspi * $duration ))
min=$(( $nspi + $buffer + $buffer ))
JOBOPT="--min $min $JOBOPT"
case "$PG" in
	-1) plus=$(( $duration + $buffer ))
	    JOBOPT="--npgplus $plus $JOBOPT";;
#		echo $JOBOPT;;
		
	1) plus=$(( $duration + $buffer ))
	   JOBOPT="--pgplus $plus $JOBOPT";;
#		echo $JOBOPT;;
	0) ;;
	*) Invalid value for playground flag;;
esac

# Tell everyone what you're doing
echo  
echo Creating files for ${runid} at `date`

# If required, create job files. Default is to create file only if one doesn't
# already exist. This can be forced to overwrite existing files.

ifopair=( "" )

first="1";
for i in $@
do
	for j in $@
	do
		if [ ! "$i" = "$j" ] && [ ! `echo $excluded | grep $i$j` ] && [ ! `echo $excluded | grep $j$i` ]; then
			echo ${ifopair[@]} | grep -q $i$j
			test1=$?
			echo ${ifopair[@]} | grep -q $j$i
			test2=$?
			if [ "$test1" -eq 1 ] && [ "$test2" -eq 1 ]; then
				if [ "$first" -eq 1 ]; then
					ifopair=( "$i$j" )
					first="0"
				else
					ifopair=( "${ifopair[@]}" "${i}${j}" )
				fi
			fi
		fi
	done
done

if [ "$JOB" == 1 ]; then
	echo "Creating job files"
#	cd ../input/jobfiles/build/
	for i in ${ifopair[@]}
	do
		if [ ! -f "${i}_${runid}.txt" ] || [ $FORCEJOB -eq 1 ]; then
			if [ $gpsEndTime -lt 900000000 ]; then
				echo "=> ../input/jobfiles/build/${i}_${runid}.txt ...done"
			
				if [ ! -f "${i}_S5.txt" ]; then
					if [ ! `echo ${i} | grep -q V1` ]; then
						if [ ! -f "${i}_vetoflags.txt" ]; then
							i1=`echo ${i} | cut -c1-2`
							i2=`echo ${i} | cut -c3-4`
							echo "-- creating ${i}_vetoflags.txt"
							cat  ../input/jobfiles/build/${i1}_vetoflags.txt ../input/jobfiles/build/${i2}_vetoflags.txt > ../input/jobfiles/build/${i}_vetoflags.txt
						fi
						echo "-- segwizarding ${i}_S5.txt"
						segwizard S5 $i `cat ../input/jobfiles/build/${i}_vetoflags.txt` > ../input/jobfiles/build/${i}_S5.txt
					else
						echo "${i}_S5.txt does not exist - please create this file using the vdb"
						exit
#					#NB this will be changed soon
					fi
				fi
				if [ "$SHELL" -eq 1 ]; then
					${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/buildJobFile.perl $JOBOPT < ../input/jobfiles/build/${i}_S5.txt > ../input/jobfiles/build/${i}_${runid}.txt
				else
					${MATAPPS_TOP}/src/searches/stochastic/CrossCorr//buildJobFile_mat.perl $JOBOPT < ../input/jobfiles/build/${i}_S5.txt > ../input/jobfiles/build/${i}_${runid}.txt
				fi
			else
				if [ ! -f "${i}_vetoflags.txt" ]; then
					i1=`echo ${i} | cut -c1-2`
					i2=`echo ${i} | cut -c3-4`
					echo "-- creating ${i}_vetoflags.txt"
					echo `cat ../input/jobfiles/build/${i1}_vetoflags.txt`,`cat ../input/jobfiles/build/${i2}_vetoflags.txt` > ../input/jobfiles/build/${i1}${i2}_vetoflags.txt
				fi
				echo "ligolw_segment_querying"
	#Choose the appropriate science mode flag for s6
	#this is still hard-coded as I've run out of params to pass to the
	#script! This may have to change as we progress with the S6 analysis
				if [ "$i1" == "V1" ]; then
					sciflag1=ITF_SCIENCEMODE;
				else
					sciflag1=DMT-SCIENCE;
				fi
				if [ "$i2" == "V1" ]; then
                                        sciflag2=ITF_SCIENCEMODE;
                                else
                                        sciflag2=DMT-SCIENCE;
                                fi
				ligolw_segment_query -t $S6_SEGMENT_SERVER -q -s $gpsStartTime -e $gpsEndTime -a ${i1}:${sciflag1},${i2}:${sciflag2} -b `cat ../input/jobfiles/build/${i}_vetoflags.txt` -o ../input/jobfiles/build/${i}_${runid}.xml
				#print to job file
				if [ "$SHELL" -eq 1 ]; then
					ligolw_print -t segment -c start_time -c end_time -d " " ../input/jobfiles/build/${i}_${runid}.xml | awk '{print "1 " $1 " " $2 " " $2-$1}' | ${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/buildJobFile.perl $JOBOPT  > ../input/jobfiles/build/${i}_${runid}.txt
				else
					ligolw_print -t segment -c start_time -c end_time -d " " ../input/jobfiles/build/${i}_${runid}.xml | awk '{print "1 " $1 " " $2 " " $2-$1}' | ${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/buildJobFile_mat.perl $JOBOPT  > ../input/jobfiles/build/${i}_${runid}.txt
				fi
			fi

			if [ ! $? == 0 ]; then
				echo "Job files not created successfully (see error message above)"
				exit 1
			fi
			echo "=> ../input/jobfiles/build/${i}_${runid}.txt ...done"
			cp ../input/jobfiles/build/${i}_${runid}.txt ../input/jobfiles/
		else
			echo "=> ../input/jobfiles/build/${i}_${runid}.txt ...file exists - skipping"
		fi
	done

#	cd ../../../command/
else
	echo "Not creating job files";
fi

# If required, create paramfiles. Default is to create file if one doesn't exist# and not to if it does. This can be forced to overwrite existing files.
if [ "$PARAM" == 1 ]; then
        echo "Running create_param_files.perl $PARAMOPT $*"
        ${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/create_param_files.perl $PARAMOPT $*

        if [ ! $? == 0 ]; then                 echo "Param files not created successfully (see error message above)"
                exit 1
        fi
else
        echo "Not creating param files";
fi


# If required, create cache files. Default is to create file only if one doesn't
# already exist. This can be forced to overwrite existing files. We are using a script that also created dag files, but we don't want to use them

if [ "$CACHE" == 1 ]; then
	echo "Creating cache files"
	for i in ${ifopair[@]}
	do
		${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/stochastic_pipe.tclsh -p ../input/paramfiles/${i}_${runid}_params.txt -j ../input/jobfiles/${i}_${runid}.txt -d ${i}_stoch_extra.dag $CACHEOPT
		#stochastic_pipe.tclsh creates a dag file that we don't want to use
		if [ ! $? == 0 ]; then
			echo "Cache files not created successfully (see error message above)"
#		#	exit 1
		fi

#		rm -f ${i}_stoch_extra.dag
	done
else
	echo Not creating cache files
fi

# If required, create dag files. Default is to create file only if one doesn't
# already exist. This can be forced to overwrite existing files.
DAGFILE="stochastic_${runid}_condor.dag"
if [ "$DAG" == 1 ]; then
	if [ ! -f "$DAGFILE" ] || [ $FORCEDAG -eq 1 ]; then
		if [ "$RENORM" == 1 ]; then
			for pair in ${ifopair[@]}; do
				modifyFilterFile="${pair}${modifyFilterSuffix}"
				if [ -z "$modifyFilter" ]; then
					modifyFilter=`cat $modifyFilterFile`
				else
					tmp=`cat $modifyFilterFile`
					modifyFilter="${modifyFilter}:${tmp}"
				fi
			done
			DAGOPT="-M $modifyFilter ${DAGOPT}"
		fi
		${MATAPPS_TOP}/src/searches/stochastic/CrossCorr/create_dag_multi.perl $DAGOPT $*
		if [ ! $? == 0 ]; then
			echo "Dag file not created successfully (see error message above)"
			exit 1
		fi

	else
		echo "Creating dag file"
		echo "=> $DAGFILE ...file exists - skipping"
	fi
else
	echo Not creating dag files
fi

# submit the dags to the cluster (unless flagged to create only).

if [ "$SUBMIT" == 1 ]; then
	condor_submit_dag $DAGFILE
else 
	echo Created files - not submitting
fi


