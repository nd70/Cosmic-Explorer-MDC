function [] = check_files_SpH()

% This code checks to see if the naive sigma files cover the same segments as 
% the SpHSet files.  The code outputs a list of jobs that need to be rerun to
% generate new naivesigma files.
clear all;

indir = '/archive/home/ethrane/S5H1L1_sph52_SID_L20/'; %files moved here-------
postfix = 'S5H1L1_sph52_SID_L20';
fileprefix = [indir postfix];
outdir = '/usr1/ethrane/';                             %files created here-----
numJobs=18837;
for ii=1:numJobs
  filename=[fileprefix '_SpH.job' num2str(ii) '.trial1.mat'];
  gpsfile=['~/S5H1L1_sph52_SID_L20/S5H1L1_sph52_SID_gps_naivesigmas.job'... 
             num2str(ii) '.trial1.dat'];
  try 
    index = load(filename);
    time = index.Sky{end}.time;
  catch 
    fprintf('%i 1 -2\n',ii);   % SpH file fails
  end
  try
    gpsData = load(gpsfile);
    time2 = gpsData(end,1);
  catch
    fprintf('%i 0 -2\n',ii);   % naivesigma file fails
  end
  try
    if time>time2
      fprintf('%i 0 -1\n',ii); % naivesigma file incomplete
    elseif time<time2
      fprinf('%i 1 -1\n',ii);  % SpH file incomplete
    else
    end
  catch
    fprintf('%i 2\n',ii);   % mystery error
  end

end % loop over ii (jobs)

fprintf('Done.\n');
return
