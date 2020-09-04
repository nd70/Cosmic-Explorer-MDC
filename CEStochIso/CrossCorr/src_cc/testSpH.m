
detectorPair='HL';
Lmax=2;
flow=40;
fhigh=1800;
deltaF=0.25;
numFreqs=(fhigh-flow)/deltaF+1;
w1w2bar=1;
w1w2squaredbar=1;
FreqIntFlag=true;

H=ones(1,numFreqs);
segmentDuration=60;
FreqIntFlag=true;
outputFilePrefix='/home/sballmer/testing/Test_A';
maxSegmentsPerFile=60;
jobNr=1;
trialNr=1;

try
  mask;
catch
  mask = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
end


vSpH=SpH(detectorPair,Lmax,numFreqs,flow,deltaF,H,mask.data,w1w2bar,w1w2squaredbar,segmentDuration,FreqIntFlag,outputFilePrefix,maxSegmentsPerFile,jobNr,trialNr);




C=zeros(1,numFreqs);
P1=ones(1,numFreqs);
P2=P1;

t=800000000+(0:60:660);
for ii=1:length(t)
  GPSstart=t(ii);
  fprintf('Calling doSpH for t=%d\n',GPSstart);
  vSpH=doSpH(vSpH,params,C,P1,P2,GPSstart);
end
vSpH=saveSpHSet(vSpH);

%N=size(glm.data,1);
%gammaS=zeros(N,N);
%XS=zeros(N,1);
%for ii=1:length(t)
%  tsidereal=t(ii);
%  gamma=calFisher(P1,P2,H,glm,tsidereal,segmentDuration,w1w2bar,w1w2squaredbar,mask.data,FreqIntFlag);
%  X=calX(C,P1,P2,H,glm,tsidereal,segmentDuration,w1w2bar,w1w2squaredbar,mask.data,FreqIntFlag);
%  gammaS=gammaS+gamma;
%  XS=XS+X;
%end
