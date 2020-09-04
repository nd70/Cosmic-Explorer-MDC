% script for testing simulateSkyMapTimeDomain and simulateSkyMapFreqDomain 

% assign some standard values to parameters for SpH search
flow = 40; % Hz
fhigh = 1800; % Hz 
deltaF = 0.25; % Hz
numFreqs = 1+floor((fhigh-flow)/deltaF);
F = flow + deltaF*(0:numFreqs-1)';

siderealtime = 0;
duration = 60; % sec
tlow = 0; % gps start time
fsample = 4096; % Hz

gammaLM_coeffsPath='/Users/joeromano/src/matapps/src/searches/stochastic/CrossCorr/';
Lmax = 5;
isSpH = true;

% get detector geometry information
det1 = getdetector('LHO');
det2 = getdetector('LLO');

% calculate gammaLMs appropriate for power spectra simulations
g1lm = calGammaLM(gammaLM_coeffsPath, 'HH', Lmax, numFreqs, flow, deltaF);
g2lm = calGammaLM(gammaLM_coeffsPath, 'LL', Lmax, numFreqs, flow, deltaF);
glm  = calGammaLM(gammaLM_coeffsPath, 'HL', Lmax, numFreqs, flow, deltaF);

% load signal power spectrum H(f) from file
Hf = load('/Users/joeromano/src/sgwb/S5/input/simfiles/HSource_1e_0.txt');
intLog = 1;

% read in a skymap
filename = '/Users/joeromano/src/sgwb/S5/input/simfiles/monopole_plmComplex_L5.txt'
filename = '/Users/joeromano/src/sgwb/S5/input/simfiles/dipole_plmComplex_L5.txt'
filename = '/Users/joeromano/src/sgwb/S5/input/simfiles/quadrupole_plmComplex_L5.txt'
filename = '/Users/joeromano/src/sgwb/S5/input/simfiles/octupole_plmComplex_L5.txt'
filename = '/Users/joeromano/src/sgwb/S5/input/simfiles/pointSource_ra6_dec45_plmComplex_L5.txt'
fileType = 2;
numBins = 1;
InjectAsSpH = true;
conversionLmax = 5;
conversionDeg = 1;
[map,coord1,coord2]=readMap(filename,fileType,numBins,InjectAsSpH,conversionLmax,conversionDeg);

% construct trivial transfer functions
transfer1 = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
transfer2 = constructFreqSeries(ones(numFreqs,1), flow, deltaF);

% calculate power spectra
[P1, P2, CP] = simulateSkyMapFreqDomain(Hf, intLog, ...
                                        flow, deltaF, numFreqs, ...    
                                        siderealtime, ...
                                        glm, g1lm, g2lm, ...
                                        det1, det2, ...
                                        isSpH, coord1, coord2, map);

%max(imag(P1))
%max(imag(P2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate time-series parameters
% NOTE: need to use different flow, deltaF, and numFreqs for gamma_LM calculation
% (the params below allow for coarser frequency resolution for the gammaLMs)
deltaF = 1/duration * 60;
N = fsample/deltaF;
% discrete positive frequencies (not including DC component and Nyquist)
if ( mod(N,2)== 0 )
  numFreqs = N/2+2;
else
  numFreqs = (N-1)/2+3;
end
flow = 1/duration;

nResample = 10;
betaParam = 5;
 
% calculate gammaLMs appropriate for time-series simulations
fprintf('calculating gamma_LMs for detector 1\n');
g1lm = calGammaLM(gammaLM_coeffsPath, 'HH', Lmax, numFreqs, flow, deltaF);
fprintf('calculating gamma_LMs for detector 2\n');
g2lm = calGammaLM(gammaLM_coeffsPath, 'LL', Lmax, numFreqs, flow, deltaF);
fprintf('calculating gamma_LMs for detectors 1 and 2\n');
glm  = calGammaLM(gammaLM_coeffsPath, 'HL', Lmax, numFreqs, flow, deltaF);

[d1, d2] = simulateSkyMapTimeDomain(duration, tlow, ...
                                    fsample, fsample, ...
                                    nResample, nResample, betaParam, betaParam, ...
                                    Hf, intLog, ...
                                    siderealtime, ...
                                    glm, g1lm, g2lm, ...
                                    det1, det2, ...
                                    isSpH, coord1, coord2, map, ...
                                    transfer1, transfer2);

% psd estimation parameters
nfft = (1/32)*fsample;
window = hann(nfft);
noverlap = nfft/2;

% estimate power spectra from time-series data
[p1,f]=pwelch(d1,window,noverlap,nfft,fsample);
[p2,f]=pwelch(d2,window,noverlap,nfft,fsample);
[cp,f]=cpsd(d1,d2,window,noverlap,nfft,fsample);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots
figure(1);
subplot(3,1,1)
plot(f,p1,'b',F,real(P1),'r');
title('Monopole');
title('Dipole');
title('Quadrupole');
title('Octupole');
title('Point source');
ylabel('P1 (1/Hz)');
legend('from time-series', 'expected', 'Location', 'SouthEast');
subplot(3,1,2)
plot(f,p2,'b',F,real(P2),'r');
ylabel('P2 (1/Hz)');
legend('from time-series', 'expected', 'Location', 'SouthEast');
subplot(3,1,3)
semilogy(f,abs(cp),'b',F,abs(CP),'r');
%ylim([1e-5 1e0]);
ylabel('CP (1/Hz)');
legend('from time-series', 'expected', 'Location', 'SouthEast');
xlabel('freq (Hz)');

fname = 'monopoleTest';
fname = 'dipoleTest';
fname = 'quadrupoleTest';
fname = 'octupoleTest';
fname = 'pointSourceTest';
print('-depsc2', fname);

