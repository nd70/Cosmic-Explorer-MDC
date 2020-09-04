function vSpH=SpH(gammaLM_coeffsPath,detectorPair,Lmax,numFreqs,flow,deltaF,H,mask,w1w2bar,w1w2squaredbar,segmentDuration,FreqIntFlag,outputFilePrefix,maxSegmentsPerFile,jobNr,trialNr)

try
  trialNr;
catch
  trialNr = 1;
end;

% constructor for the Spherical Harmonics class
%

vSpH.gammaLM_coeffsPath = gammaLM_coeffsPath;
vSpH.detectorPair=detectorPair;
vSpH.Lmax=Lmax;
vSpH.numFreqs=numFreqs;
vSpH.flow=flow;
vSpH.deltaF=deltaF;
vSpH.H=H;
vSpH.mask=mask;
vSpH.w1w2bar=w1w2bar;
vSpH.w1w2squaredbar=w1w2squaredbar;
vSpH.segmentDuration=segmentDuration;
vSpH.FreqIntFlag=FreqIntFlag;

vSpH.prefix   =[outputFilePrefix,'_SpH'];
vSpH.suffix   =['.job' num2str(jobNr) '.trial' num2str(trialNr) '.mat'];
vSpH.maxSegmentsPerFile=maxSegmentsPerFile;
vSpH.setPrefix=[outputFilePrefix,'_SpHSet'];
vSpH.setNumber=1;
vSpH.setOffset=0;
vSpH.currentFilename=[vSpH.setPrefix,num2str(vSpH.setNumber),vSpH.suffix];

vSpH.out.counter=0;
vSpH.out.data={};
vSpH.out.meta={};

vSpH.glm=calGammaLM(vSpH.gammaLM_coeffsPath,vSpH.detectorPair,vSpH.Lmax,vSpH.numFreqs,vSpH.flow,vSpH.deltaF);
vSpH.g1lm=calGammaLM(vSpH.gammaLM_coeffsPath,[vSpH.detectorPair(1),vSpH.detectorPair(1)],vSpH.Lmax,vSpH.numFreqs,vSpH.flow,vSpH.deltaF);
vSpH.g2lm=calGammaLM(vSpH.gammaLM_coeffsPath,[vSpH.detectorPair(2),vSpH.detectorPair(2)],vSpH.Lmax,vSpH.numFreqs,vSpH.flow,vSpH.deltaF);

vSpH = class(vSpH,'SpH');

