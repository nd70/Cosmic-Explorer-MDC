% script for testing simulateDetectorNoise

% assign some standard values to parameters for SpH search
flow = 40; % Hz
fhigh = 1800; % Hz 
deltaF = 0.25; % Hz
numFreqs = 1+floor((fhigh-flow)/deltaF);
F = flow + deltaF*(0:numFreqs-1)';

duration = 600; % sec
tlow = 0; % gps start time
fsample = 4096; % Hz
nResample = 10;
betaParam = 5;

% load noise power spectra P1(f), P2(f) from file
%P1 = load('/home/sballmer/stochastic/sgwb/S5/input/simfiles/LIGOsrdPSD.txt');
%P2 = load('/home/sballmer/stochastic/sgwb/S5/input/simfiles/LIGOsrdPSD.txt');
P1 = load('/Users/joeromano/src/sgwb/S5/input/simfiles/LIGOsrdPSD.txt');
P2 = load('/Users/joeromano/src/sgwb/S5/input/simfiles/LIGOsrdPSD.txt');
intLog = 1;

% construct trivial transfer functions
transfer1 = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
transfer2 = constructFreqSeries(ones(numFreqs,1), flow, deltaF);

if 0
[d1, d2] = simulateDetectorNoise(duration, tlow, ...
                                 fsample, fsample, ...
                                 nResample, nResample, ...
                                 betaParam, betaParam, ...
                                 P1, P2, intLog, ...
                                 transfer1, transfer2);
else
  simDetectorNoisePowerSpec1=P1;
  simDetectorNoisePowerSpec2=P2;
  params.resampleRate1=fsample;
  params.resampleRate2=fsample;
  params.nResample1=nResample;
  params.nResample2=nResample;
  params.betaParam1=betaParam;
  params.betaParam2=betaParam;
  params.simulatedSkyMapInterpolateLogarithmic=true;
  params.flow=flow;
  params.deltaF=deltaF;
  params.numFreqs=numFreqs;
  params.ASQchannel1='';
  params.alphaBetaFile1='none';
  params.calCavGainFile1='none';
  params.calResponseFile1='none';
  params.cal1.t=0; params.cal1.f=0; params.cal1.R0=0; params.cal1.C0=0; params.cal1.alpha=0; params.cal1.gamma=0;
  params.ASQchannel2='';
  params.alphaBetaFile2='none';
  params.calCavGainFile2='none';
  params.calResponseFile2='none';
  params.cal2=params.cal1;
  params.simulatedSkyMapBufferDepth=6000;
  params.simulatedSkyMapHalfRefillLength=32;
  initDetectorNoiseData(params.resampleRate1,params.resampleRate2,...
                 params.nResample1,params.nResample2,...
                 params.betaParam1,params.betaParam2,...
                 simDetectorNoisePowerSpec1,simDetectorNoisePowerSpec2,...
                 params.simulatedSkyMapInterpolateLogarithmic,...
		 params.flow,params.deltaF,params.numFreqs,...
                 params.cal1.t,params.cal1.f,params.cal1.R0,params.cal1.C0,params.cal1.alpha,params.cal1.gamma,...
		 params.ASQchannel1,params.alphaBetaFile1,params.calCavGainFile1,params.calResponseFile1,...
		 params.cal2.t,params.cal2.f,params.cal2.R0,params.cal2.C0,params.cal2.alpha,params.cal2.gamma,...
		 params.ASQchannel2,params.alphaBetaFile2,params.calCavGainFile2,params.calResponseFile2,...
		 params.simulatedSkyMapBufferDepth,params.simulatedSkyMapHalfRefillLength);
GPSstart=900000000;
dur=duration;
[d1, d2, badResponse] = getDetectorNoiseData(GPSstart,dur);
end				 

% psd estimation parameters
Nt = length(d1);
nfft = Nt/160;
window = hann(nfft);
noverlap = nfft/2;

% estimate power spectra from time-series data
[p1,f1]=pwelch(d1,window,noverlap,nfft,fsample);
[p2,f2]=pwelch(d2,window,noverlap,nfft,fsample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
figure(1)
plot(d1)
xlabel('sample');
ylabel('time series data');
title('Detector 1')

figure(2)
plot(d2)
xlabel('sample');
ylabel('time series data');
title('Detector 2')

figure(3);
loglog(f1,p1,'b',P1(:,1),P1(:,2),'r');
xlim([F(1) F(end)]);
ylabel('P (strain^2/Hz)');
legend('from time-series', 'expected', 'Location', 'SouthEast');
title('Detector 1')

figure(4);
loglog(f2,p2,'b',P2(:,1),P2(:,2),'r');
xlim([F(1) F(end)]);
ylabel('P (strain^2/Hz)');
legend('from time-series', 'expected', 'Location', 'SouthEast');
title('Detector 2')

% coherence parameters
noverlap = Nt/16;
window = hann(nfft);
noverlap = nfft/2;
[cp,f]=mscohere(d1,d2,window,noverlap,nfft,fsample);

figure(5);
loglog(f,abs(cp));
xlim([F(1) F(end)]);
ylabel('coherence');
xlabel('freq (Hz)');

