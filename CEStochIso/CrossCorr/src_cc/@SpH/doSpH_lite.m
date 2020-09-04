function vSpH=doSpH_lite(vSpH,params,C,P1,P2,GPSstart)

% This function does the logistics of the sperical harmonics calculation
%
%

% NEVER SAVE SPH_SET FILES!!!
%if vSpH.out.counter-vSpH.setOffset >= vSpH.maxSegmentsPerFile
%  vSpH=saveSpHSet(vSpH,false); %SpHSet is full, update vSpH object 
%end

% Calculate sidereal time as it will be used for calculations involving gamma.
GPSmid=GPSstart+vSpH.segmentDuration/2;
tsidereal=GPStoGreenwichMeanSiderealTime(GPSmid);

ii=vSpH.out.counter+1; %meta data index
jj=ii-vSpH.setOffset;  %vSpH index

% Fill vSpH object with data such as X, Fisher, coherence, etc..
vSpH.out.data{jj}.time=GPSstart;
vSpH.out.data{jj}.midSegSidereal=tsidereal;
vSpH.out.data{jj}.X      = calX   (C,P1,P2,vSpH.H,vSpH.glm,tsidereal,vSpH.segmentDuration,vSpH.w1w2bar,vSpH.w1w2squaredbar,vSpH.mask,vSpH.FreqIntFlag);
vSpH.out.data{jj}.Fisher = calFisher(P1,P2,vSpH.H,vSpH.glm,tsidereal,vSpH.segmentDuration,vSpH.w1w2bar,vSpH.w1w2squaredbar,vSpH.mask,vSpH.FreqIntFlag);
[coh_f,coh] = calCoherence(C,P1,P2,vSpH.deltaF,vSpH.segmentDuration,vSpH.w1w2bar,vSpH.w1w2squaredbar,vSpH.mask); 
vSpH.out.data{jj}.coh = coh;

% only write coh_f to files if desired (coh used for model selection)
if params.writeCohFToFiles
  vSpH.out.data{jj}.coh_f = coh_f;
end

% calculate Psi and write to files only if simulated sky map in SpH basis
% Psi is a diagnostic tool used to investigate small signal approx..
if params.doSimulatedSkyMap && params.simulatedSkyMapInjectAsSpH
  vSpH.out.data{jj}.Psi = calPsi(P1,P2,vSpH.H,vSpH.glm,tsidereal,vSpH.segmentDuration,vSpH.w1w2bar,vSpH.w1w2squaredbar,vSpH.mask,vSpH.FreqIntFlag);
end

% Update metadata
vSpH.out.meta{ii}.time=GPSstart;
vSpH.out.meta{ii}.midSegSidereal=tsidereal;
vSpH.out.meta{ii}.filename=vSpH.currentFilename;
vSpH.out.meta{ii}.segmentOffset=vSpH.setOffset;
vSpH.out.counter=ii;

% Notes: use array of struct instead of cell of struct - conversion see below
% conversion between and new file format:
% B=p.Sky;
% A=cell2mat(B);
% B=mat2cell(A,1,ones(1,length(A)));
