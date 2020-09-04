function params = addSTAMPparams(params)
% Used when running preproc directly from stochmap.
% Used to set default parameters if they aren't already set.
% T. Prestegard
% Last updated 4/25/2011

% Structure of default STAMP parameters.
default.G=6.673e-11;
default.yMapScale=1e-43;
default.c=299792458;
default.doMC=false;
default.doRadon=false;
default.doPolar=false;
default.purePlus=false;
default.pureCross=false;
default.doBoxSearch=false;
default.doLH=false;
default.doRadiometer=false;
default.saveMat=false;
default.savePlots=true;
default.debug=false;
default.doStampFreqMask=1;
default.StampFreqsToRemove = [];
default.StampnBinsToRemove = [];
default.alphaExp=[];
default.kludge=1;
default.FMapScale=1e3;
params.ErgsPerSqCM=1000;
params.fluence=0;
default.fft1.dataWindow=-1;
default.fft2.dataWindow=-1;
default.doRadonReconstruction=false;
default.doClusterSearch=false;
default.fixAntennaFactors=false;
default.phaseScramble=false;
default.glitchCut=false;
default.bknd_study = false;
default.bknd_study_dur = 10;
default.Autopower = false;
default.pixelScramble = false;
default.NoTimeDelay = false;
default.recordCuts = false;
default.coarsegrained = false;
default.stamp_pem = false;
default.DQmatfile = 'input';
default.skypatch = false;
default.fastring = false;
default.outputfilename = 'map';
default.discontinuity = false;
default.powerinj = false;
default.inj_trials = 1;
default.alpha_n=1;
default.bkndstudy=false;
default.twojoblimit=true;
default.doCoincidentCut=false;
default.loudPixel=false;
default.crunch_map=false;
default.doBurstegard=false;

names=fieldnames(default);

for i=1:length(names)
  if ~isfield(params,names{i})
    params=setfield(params,names{i},getfield(default,names{i}));
  end
end

% fft substruct
subnames1=fieldnames(default.fft1); subnames2=fieldnames(default.fft2);
for i=1:length(subnames1)
  if ~isfield(params.fft1,subnames1{i})
    params.fft1=setfield(params.fft1,subnames1{i},getfield(default.fft1,subnames1{i}));
  end
end
for i=1:length(subnames2)
  if ~isfield(params.fft2,subnames2{i})
    params.fft2=setfield(params.fft2,subnames2{i},getfield(default.fft2,subnames2{i}));
  end
end

return
