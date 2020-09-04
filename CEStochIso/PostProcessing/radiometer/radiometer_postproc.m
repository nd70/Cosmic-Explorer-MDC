function radiometer_postproc(paramsFile,post_proc_paramsFile,jobsFile,epoch,jobNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%paramsFile --------------> paramsFile used in radiometer analysis.
%post_proc_paramsFile ====> post processing params file containing save location, data quality cut information, etc.
%epoch -------------------> analysis epoch that you wish to run this over.
%jobNumber ---------------> ceil(numJobs/250). Used for running this on condor.
%This is a post-processing script for the radiometer code
%It finds the bad GPS times and sets up windowing information, then calls combineResultsFromMultipleJobsRM.m
%
%
% Written by stefan ballmer, adapted by patrick meyers
% meyers@physics.umn.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
% for condor purposes
%
%%%%%%%%
if isstr(jobNumber)
jobNumber        = str2num(jobNumber);
end

pproc_params = readParamsFromFile_radiometer_post_proc(post_proc_paramsFile);
pproc_params = pprocParamsCheck(pproc_params);

if strcmp(epoch,'SumLocalRad')
  if pproc_params.doNarrowbandRadiometer
    SumLocalRad(post_proc_paramsFile);
  else
    SumLocalRad_skymaps(post_proc_paramsFile);
  end
  return;
end

% load in params used in stochastic analysis
%
params           = readParamsFromFile(paramsFile);

% load in params to be used for post processing analysis

% define output file prefix from stochastic analysis params file
outputFilePrefix = params.outputFilePrefix;

% define segment duration
segmentDuration  = params.segmentDuration;

% H1L1, for example is used consistently throughout the post processing code
ifopair = [params.ifo1 params.ifo2];



% for saving purposes and purposes of the cut. Currently the file prefix (absent the preceeding directory structure) has the form: (ifopair)_S6_blnd_(epoch)_(flow)_(fhigh)
%temp  = regexp(outputFilePrefix, '/', 'split');
%temp2 = regexp(temp{end},'_','split');

fhigh = params.fhigh;
flow  = params.flow;
save_loc = pproc_params.post_processing_outputFilePrefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileSuffix='.trial1.mat';
fprintf('loading files...\n');
cjobs=[1:pproc_params.jobsToCombine]+(jobNumber-1)*pproc_params.jobsToCombine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stationarity cut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badGPSTimes = [];

% initiate Bad GPS Times vector
if exist(pproc_params.badGPSTimesFile)
    [val st et dur] = textread(jobsFile,'%d %d %d %d');
    % get start time / end time for cjobs
    st = st(cjobs(1));
    end_index = min(length(et),cjobs(end));
    et = et(end_index);
    [badGPSTimes_temp ifopairs flows fhighs] = textread(pproc_params.badGPSTimesFile,'%f %s %f %f');
    for ii = 1:length(badGPSTimes_temp)
        if and(flows(ii) == flow, fhighs(ii) == fhigh)
            if strcmp(ifopairs{ii}, ifopair)
                if and(badGPSTimes_temp(ii) > st, badGPSTimes_temp(ii) < et)
                  badGPSTimes = [badGPSTimes;badGPSTimes_temp(ii)];
                end
            end 
        end 
    end 
end
totSeg=0;
lenp=0;
fileGood=1;
% do delta sigma cut on the fly
if pproc_params.doDeltaSigmaCut
    dsigvetothresh = pproc_params.deltaSigmaThreshold; % delta sigma cut threshold. 
    index=0;
    for ii=cjobs
        index=index+1;
        %find badGPStimes for this specific job
        fileGood=1;
    try
        file = [outputFilePrefix '_naivesigmas.job' num2str(ii) '.trial1.dat'];
        p=load([outputFilePrefix '_naivesigmas.job' num2str(ii) '.trial1.dat']);
    catch p=[]; end;%trycatch
    lenp=length(p);
    if lenp==0
      fileGood=0;
    end;%lenp
    if fileGood
      %stationarity cut determines bad GPS times, see ballmer thesis
      bad = find(~(abs(log(p(:,2)./p(:,3))) < log(dsigvetothresh)));
      badGPSTimes=[badGPSTimes;p(bad,1)];
      totSeg=totSeg+lenp;
    else
      fprintf('Could not load file %s %d\n', 'naive sigma job', ii);
    end % fileGood
  end % for ii=cjobs
end % if do delta sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up windowing information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~(params.hannDuration1==segmentDuration)|~(params.hannDuration2==segmentDuration)
  fprintf('WARNING Hann Window duration is not equal to segment duration\n');
end% warning
hannDuration1=params.hannDuration1;
hannDuration2=params.hannDuration2;
if params.doOverlap
  fprintf('Overlapping segments is turned on. We use pure Hann Window\n');
  hannDuration1=segmentDuration;
  hannDuration2=segmentDuration;
end%do overlap hann hann window settings

numPoints1=segmentDuration*params.resampleRate1;
numPoints2=segmentDuration*params.resampleRate2;
window1=tukeywin(numPoints1,hannDuration1/segmentDuration);
window2=tukeywin(numPoints2,hannDuration2/segmentDuration);

fprintf('combining data with output file prefix = %s\n number of jobs =%d\n segment duration =%d\n number of bad GPS times=%d\n',outputFilePrefix,length(cjobs),segmentDuration,length(badGPSTimes));

    [sky,numSegments]=combineResultsFromMultipleJobsRM_2(...
                      outputFilePrefix, ...
                      fileSuffix,...
                      cjobs,...
                      segmentDuration,...
                      badGPSTimes,...
                      params.doOverlap,...
                      window1,window2);

    fprintf('number of segments = %d\n bad GPS Times Percent: %4.2f\n',numSegments,(length(badGPSTimes)-1)/(length(badGPSTimes)+numSegments - 1)*100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS
% results saved in format:
% (post_proc_outputFilePrefix)_(epoch)_combined_jobs_(jobNumber).trial1.mat
% where, again, jobNumber specifies which 250 jobs have been combined for the stochastic jobs
% because of the way the condor dags are set up it's easier to put epoch as an input for this than to
% specify it as a parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([save_loc '_' ifopair '_' epoch '_' num2str(flow) '_' num2str(fhigh) '_combined_jobs_' num2str(jobNumber)  fileSuffix], 'sky', ...
     'numSegments', 'badGPSTimes','paramsFile','jobsFile','flow','fhigh','pproc_params');


