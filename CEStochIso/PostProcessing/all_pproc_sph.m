function all_pproc_sph(str, lmax, outdir, matdir, matstr, badgpsstr, cuts, njobs)
%function   all_pproc_sph(str, lmax, outdir, matdir, badgpsstr, cuts)
%
% combines outpt jobs for sph analysis
%
% str       : string specifying run parameters in output filename
% lmax      : resolution cut-off used during analysis
% outdir    : directory where output data is located
% matdir    : directory where .mat file should be stored
% matstr    : string for naming output .mat file            
% badgpsstr : file string for bad gps time to remove
% cuts      : set this >0 if you want both combined naive sigma and badgps cuts
%
% Author: E Thrane with modifications by L. Sammut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmax = strassign(lmax);

 % check if output directory is specified, defaults to working/output directory
 if ~exist('outdir','var')
      outdir = './output';
 end
 % check if output directory is specified, defaults to working directory
 if ~exist('matdir','var')
      matdir = '.';
 end

 try
   cuts;
 catch
   cuts = 0;
 end;

% check if badgpsstr is specified, otherwise use naive sigmas
   if exist('badgpsstr','var') && cuts>0 ;
     badgpsswitch = 1;
     try
       badgps=load(badgpsstr);, % load bad gps times
     catch
       badgpsswitch = 0;, 
     end
   elseif exist('badgpsstr','var')
     badgpsswitch = 2;
     try
       badgps = load(badgpsstr);, % load bad gps times
     catch
       badgpsswitch = 0;,
     end 
   else
     badgpsswitch = 0 ;
   end

% for overlapping Hann-windowed segments
epsilon=3/70;

% number of parameters
N=(lmax+1)^2;

% intitalize X, Fisher, etc..
     X=zeros(N,1); Fisher=zeros(N,N); bad_segs=0; bad_jobs=0; nsegs=0; 

% loop over jobs combine SkySet files for each job to create a single SKY array
files = dir([outdir '/' str '_SpHSet*']);

% check for number of jobs input
try
  njobs;
catch
  njobs = 2*length(files);
end;

for ii=1:njobs; %18837
  try % this is to catch any unaccounted for fails
    SKY = [];
    for jj = 1:20  % number of possible SpHSet files
      matfile = ...
        [outdir '/' str '_SpHSet' num2str(jj) '.job' num2str(ii) '.trial1.mat'];
      try
        g = load(matfile);
        SKY = [SKY g.Sky];
      catch
        % do nothing
      end
    end

    % initialize fisher1 and fisherN for endpoint effects
    fisher1=zeros(N,N); fisherN=zeros(N,N);

    % read in naivesigma data
    naifile = [outdir '/' str '_naivesigmas.job' num2str(ii) '.trial1.dat'];
    naidat = load(naifile);
    [nrows ~] = size(naidat);

    if nrows~=length(SKY)
      fprintf('mismatch with job %i\n', ii);
      bad_jobs = bad_jobs + 1;
    elseif length(SKY)==0
      fprintf('Job %i does not contain any data \n', ii);
      bad_jobs = bad_jobs + 1;
    else % not a bad job so proceed with calculation
      lastsegmentwasbad=1;
      for kk=1:length(SKY)
        % stationarity cut using naive sigmas or bad gps times switch
	switch badgpsswitch
          case 0
               % use naive sigmas threshold of 20%
	       cutif = naidat(kk,2)/naidat(kk,3)<0.8 | naidat(kk,2)/naidat(kk,3)>1.2;
          case 1 
               % use naive sigmas threshold of 20% AND bad GPS times
	       cutif = naidat(kk,2)/naidat(kk,3)<0.8 | naidat(kk,2)/naidat(kk,3)>1.2 | ismember(naidat(kk,1),badgps) > 0;
          case 2
               % compare bad gps time instead of naive sigmas
               cutif = ismember(naidat(kk,1),badgps) > 0; % check whether segment time is in badgps times
        end
        if cutif
	   bad_segs = bad_segs +1; % bad segment
           lastsegmentwasbad=1 ;
        else
          X = X + SKY{kk}.X;
          Fisher = Fisher + SKY{kk}.Fisher;
          nsegs=nsegs+1;
          if lastsegmentwasbad % endpoint effect
	    Fisher = Fisher+((epsilon/2)/(1-2*epsilon))*(fisher1+fisherN);
            fisher1 = SKY{kk}.Fisher;
          end
	  fisherN = SKY{kk}.Fisher;
          lastsegmentwasbad=0;
        end
      end
    end
  catch
    fprintf('other bad job: %i\n',ii);
  end
end

fprintf('bad segs = %i\n',bad_segs);
fprintf('bad jobs = %i\n',bad_jobs);
fprintf('good segs = %i\n',nsegs);

% epsilon scaling
X = (1-2*epsilon)*X;
Fisher = (1-2*epsilon)*Fisher + (epsilon/2)*(fisher1+fisherN);

save([matdir '/all_pproc_sph_' matstr '.mat'], 'X', ...
  'Fisher', 'bad_segs', 'bad_jobs', 'nsegs');
return
