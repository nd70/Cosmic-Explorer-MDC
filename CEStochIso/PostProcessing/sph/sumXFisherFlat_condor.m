function [] = sumXFisher_S5(ii)

% This script replaces sumXFisher.  It combines the data and calculates C1(t)
% By Eric Thrane.  Modified Oct 13 to skip over segments failing sigma cut.
% !!!!!!!!!!!   Must run bad.m first to generate bad_segs.txt   !!!!!!!!!!!!!

LMAX=20;

% List of segments failing sigma cut------------------------------------------
sigmacut=load('/archive/home/ethrane/SID/bad_segs.txt');
sigmacut_starts=sigmacut(:,3);
% Times identified as bad times for S5----------------------------------------
badGPS=load('/archive/home/ethrane/sgwb/S5/output/S5H1L1_badGPSsegments.dat');
bad_start=badGPS(:,1);  bad_end=badGPS(:,2);
% Jobs cut because naivesigma file could not load-----------------------------
jobcut=load('/archive/home/ethrane/SID/bad_jobs.txt');
% Example of a bad job: job 6 runs w/o errors on stochastic, but there are no
% segments in the output file.

t=1;
epsilon = 3/70;   %overlap factor from Hann windows
L0=20;            %lmax
%cet indir = '/archive/home/ethrane/S5H1L1_sph52_SID_L20/'; %files moved here--
indir = '/archive/home/ethrane/S5H1L1_sph52_flat_L20/'; %files moved here-----
%cet postfix = 'S5H1L1_sph52_SID_L20';
postfix = 'S5H1L1_sph52_SID_flatStrain';
fileprefix = [indir postfix];
outdir = '/usr1/ethrane/';                             %files created here-----
x_opt = zeros((L0+1)^2,1);
fisher_opt = zeros((L0+1)^2);
segtot=0;
coh=0;

  filename=[fileprefix '_SpH.job' num2str(ii) '.trial1.mat'];
  try 
    index = load(filename);
    Nsegs = length(index.Sky);

    % Initialize variables
    loadeddatafile='';
    dummy_file=index.Sky{1}.filename;
    %filename points to output dir, but the files have been moved.
    dummy_file=regexprep(dummy_file,outdir,indir);
    dummy=load(dummy_file);
    x=zeros(size(dummy.Sky{1}.X)); % initialize x and fisher for this job------
    fisher=zeros(size(dummy.Sky{1}.Fisher)); %---------------------------------

    for kk=1:Nsegs % loop over segments in the job-----------------------------
      datafilename = index.Sky{kk}.filename;
      % filename points to output dir, but the files have been moved.
      datafilename=regexprep(datafilename,outdir,indir);
      if strcmp(datafilename,loadeddatafile);
        jj=jj+1; % an index to keep track of where we are in a file
      else
        loadeddatafile=datafilename;
        datafile = load (datafilename);
        jj=1;
      end;

      % Oct 13: do not add any segments failing cuts---------------------------
      if any(sigmacut_starts==datafile.Sky{jj}.time)
      else % The data do not fail any cuts-------------------------------------
	coh = coh + datafile.Sky{jj}.coh;
        x = x + datafile.Sky{jj}.X;
        fisher = fisher + datafile.Sky{jj}.Fisher;
      end
    end % loop over segments kk in the job-------------------------------------
    fisher1 = datafile.Sky{1}.Fisher;
    fisherN = datafile.Sky{end}.Fisher;
    X = (1-2*epsilon)*x;
    Fisher = (1-2*epsilon)*fisher + (epsilon/2)*(fisher1+fisherN);
    time = datafile.Sky{end}.time;
    segtot=segtot+Nsegs;

  catch
    X = zeros((L0+1)^2,1);
    Fisher = zeros((L0+1)^2);
  end

  x_opt=x_opt+X;
  fisher_opt=fisher_opt+Fisher;

%cet  save(['/usr1/ethrane/sumXFisher_job' num2str(ii) '.mat'], ...
%cet    't','coh','datafilename','x_opt','fisher_opt', 'segtot');
  save(['/usr1/ethrane/sumXFisherFlat_job' num2str(ii) '.mat'], ...
    't','coh','datafilename','x_opt','fisher_opt', 'segtot');

return
