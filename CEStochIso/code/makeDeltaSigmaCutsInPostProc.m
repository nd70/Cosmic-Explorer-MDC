function makeDeltaSigmaCutsInPostProc(jobfile, stochOutputPrefix, dsc, doOverlap, segmentDuration, deltaF, outname)
%
% This script goes through the *_naivesigmas* output files for all of the 
% jobs processed by stochastic.m, calculates a user-specified delta sigma 
% cut (taking into account bias factor), and writes to file
% sGc 02/15/16
%
% Input
% - jobfile = jobfile used for particular run of stochastic.m
% - stochOutputPrefix = prefix for output of stochastic.m
% - dsc = desired delta sigma cut
% - doOverlap = 0 for no stochastic.m overlap of time segments (non-standard),
%   1 for overlap
% - segmentDuration = segment duration used for stochastic.m run
% - deltaF = frequency bin size used in stochastic.m
% - outname = name of output bad GPS times based on delta sigma cut
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%USER-SPECIFIED INPUT%%%
%Input
%jobfile='/home/sgwynne.crowder/sgwb/O1_zerolag/C02/input/jobfiles/JOB-FILE-1126623617-1136649617-master_C02-events-removed.dat';
%stochOutputPrefix='/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/a3/a3';

%Desired delta sigma cut
%dsc=0.2;

%Info for calculating bias factor for naive and sliding sigmas (info about bias factor on ilog & Stefan Ballmer's thesis)
%doOverlap=1;
%segmentDuration=192;
%deltaF=1/32;
if doOverlap
  segs=segmentDuration*deltaF*2-1;
else
  segs=segmentDuration*deltaF;
end
nn=2*9/11*segs; %9/11 for Welch factor
bf_ss=nn/(nn-1); %sliding sigma bias factor
nn=9/11*segs;
bf_ns=nn/(nn-1); %naive sigma bias factor
%bf_ss=1; %trying getting rid of bias factor
%bf_ns=1; %trying getting rid of bias factor

%Output
%outname='/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/a0_postproc/badGPSTimes_a3.txt';
%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters
jobs=load(jobfile);
num_jobs=size(jobs,1);

%Process
fid=fopen(outname,'w+');
for ii=1:num_jobs %step through jobs
  sigmasfilename=[stochOutputPrefix '_naivesigmas.job' num2str(ii) '.trial1.dat'];
  if exist(sigmasfilename)
    sigmas=load(sigmasfilename,'-ascii'); %column 1 = start sec, column 2 = naive sigma, column 3 = sliding sigma
    [m,n]=size(sigmas);
    if m~=0
      for jj=1:m %write times exceeding delta sigma cut to file
	if abs(sigmas(jj,3)*bf_ss-sigmas(jj,2)*bf_ns)/(sigmas(jj,3)*bf_ss)>dsc %apply bias factor
	  %fprintf(fid,'%i  %e  %e\n',sigmas(jj,1),sigmas(jj,2),sigmas(jj,3));
	  fprintf(fid,'%i\n',sigmas(jj,1));
        end
      end
    end
  end
end %jobs
fclose(fid)

return;
