function master_cache()
% function master_cache()
% E. Thrane on August 12, 2010 ------------------------------------------------
% This script combines all the S5 cachefiles for H1 and L1 into a single struct
% called cache which is saved to cache.mat.
%

% define paths to cachefiles
%gpsTimesPath1 = '/archive/home/shivaraj/sgwb/S5/input/cachefiles/caltech/H1L1/S5H1L1_full_run/H1/';
%gpsTimesPath2 = '/archive/home/shivaraj/sgwb/S5/input/cachefiles/caltech/H1L1/S5H1L1_full_run/L1/';
%frameCachePath1 = '/archive/home/shivaraj/sgwb/S5/input/cachefiles/caltech/H1L1/S5H1L1_full_run/H1/';
%frameCachePath2 = '/archive/home/shivaraj/sgwb/S5/input/cachefiles/caltech/H1L1/S5H1L1_full_run/L1/';
gpsTimesPath1 = '/archive/home/ethrane/cache/S5H1L1/';
gpsTimesPath2 = '/archive/home/ethrane/cache/S5H1L1/';
frameCachePath1 = '/archive/home/ethrane/cache/S5H1L1/';
frameCachePath2 = '/archive/home/ethrane/cache/S5H1L1/';


NJobs = 18837;    % total number of jobs to combine
site = ['H' 'L']; % sites for ifos to cross-correlate

% initialize cache struct
cache.frames1 = '';
cache.frames2 = '';
cache.gps1 = [];
cache.gps2 = [];

NFail = 0;        % keep track of number of failed jobs
for ii=1:NJobs
  try
    % load cachefiles
    frames1 = textread([frameCachePath1 'frameFiles' site(1) '.' num2str(ii) ...
					'.txt'],'%s');
    frames2 = textread([frameCachePath2 'frameFiles' site(2) '.' num2str(ii) ...
					'.txt'],'%s');
    gps1 = load([gpsTimesPath1 'gpsTimes' site(1) '.' num2str(ii) '.txt']);
    gps2 = load([gpsTimesPath2 'gpsTimes' site(2) '.' num2str(ii) '.txt']);

    % fill cache struct
    cache.frames1 = [cache.frames1 ; frames1];
    cache.frames2 = [cache.frames2 ; frames2];
    cache.gps1 = [cache.gps1 ; gps1];
    cache.gps2 = [cache.gps2 ; gps2];
  catch
    fprintf('Error with job %i.\n',ii);
    NFail = NFail + 1;
  end
end

% save metadata
cache.site = site;
cache.path1A = gpsTimesPath1;
cache.path2A = gpsTimesPath2;
cache.path1B = frameCachePath1;
cache.path2B = frameCachePath2;
cache.NJobs = NJobs;
cache.NFail = NFail;

%------------------------------------------------------------------------------
% clean-up: now remove duplicate gps times so that the frame2 and gps2 arrays
% contain only unique gps times
%------------------------------------------------------------------------------

last_gps2 = 0;
jj = 0;  % count duplicates
kk = 1;  % count unique frames
for ii=1:length(cache.gps2)
  if cache.gps2(ii)==last_gps2
    jj = jj + 1;   % remove duplicate frames
  else
    gps2_new(kk) = cache.gps2(ii);
    frames2_new(kk) = cache.frames2(ii);
    kk = kk + 1;
    last_gps2 = cache.gps2(ii);
    if mod(kk,1000) == 0
      fprintf('%i\n',kk);   % update progress
    end
  end  
end

cache.gps2 = gps2_new';
cache.frames2 = frames2_new';
cache.nduplicates = jj;

% save to mat file
save('./cache.mat','cache');

return

