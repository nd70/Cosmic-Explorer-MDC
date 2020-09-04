function vSpH=callSpH(vSpH,params,rbartilde1,rbartilde2,calPSD1_avg,calPSD2_avg,GPSstart)

% This prepares C, P1 and P2 and calls doSpH
%
% $Id:

global SKY_MAP_MEMORY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alot of this code was taken from (or modelled after) calCrossCorr

% extract frequency series metadata 
[data1, flow1, deltaF1, symmetry1] = extractFreqSeries(calPSD1_avg);
numFreqs1 = length(data1);

[data2, flow2, deltaF2, symmetry1] = extractFreqSeries(calPSD2_avg);
numFreqs2 = length(data2);

% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that flow, deltaF, numFreqs agree
if flow1~=flow2
  error('flow mismatch');
end 
flow = flow1;

if deltaF1~=deltaF2
  error('deltaF mismatch');
end 
deltaF = deltaF1;

if numFreqs1~=numFreqs2
  error('numFreq mismatch');
end 
numFreqs = numFreqs1;

% check that flow is >= 0
if ( flow < 0 )
  error('flow < 0');
end;

% check that frequencies actually overlap
if ( rbartilde1.flow > ...
     rbartilde2.flow + rbartilde2.deltaF * (length(rbartilde2.data)-1) ...
     | rbartilde2.flow > ...
     rbartilde1.flow + rbartilde1.deltaF * (length(rbartilde1.data)-1) )
  error('frequency ranges do not overlap');
end;

% check that the rbartildes have the same deltaF
if ( rbartilde1.deltaF ~= rbartilde2.deltaF )
  error('input data deltaF mismatch');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trim rbartildes so they refer to the same frequencies

offset = (rbartilde2.flow-rbartilde1.flow)*(1/rbartilde1.deltaF);

% check that offset in flows is an integer number of freq bins
if ( offset ~= floor(offset) )
  error('input data flow mismatch');
end;
if offset > 0
  rbartilde1.flow = rbartilde2.flow;
  rbartilde1.data = rbartilde1.data((1+offset):end);
elseif offset < 0
  rbartilde2.flow = rbartilde1.flow;
  rbartilde2.data = rbartilde2.data((1-offset):end);
end;

if length(rbartilde1.data) > length(rbartilde2.data)
  rbartilde1.data = rbartilde1.data(1:length(rbartilde2.data));
elseif length(rbartilde2.data) > length(rbartilde1.data)
  rbartilde2.data = rbartilde2.data(1:length(rbartilde1.data));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form the cross-power
data = conj(rbartilde1.data) .* rbartilde2.data;
rr   = constructFreqSeries(data, rbartilde1.flow, rbartilde1.deltaF);

% coarse grain the product to agree with power spectrum density
% calculate number of fequencies in power spectrum density 
% NOTE: index1, index2, frac1, frac2 are not needed for our analysis
[rrCG, index1, index2, frac1, frac2] = ... 
             coarseGrain(rr, params.flow, params.deltaF, params.numFreqs);
C    = 2* rrCG.data ./(vSpH.w1w2bar * params.segmentDuration);
% C corresponds to C_ft; holds normalized, coarse-grained cross-power, one-sided
P1   = calPSD1_avg.data; % power spectral density, interferometer 1
P2   = calPSD2_avg.data; % power spectral density, interferometer 2

% inject into power spectra and cross-power if desired
if (params.doSimulatedSkyMap==true) & (params.simulatedSkyMapInjectTimeDomain==false)
  [simP1, simP2, simC] = simulateSkyMapFreqDomain(SKY_MAP_MEMORY.Hf,...
                                                  SKY_MAP_MEMORY.intLog,...
	                                          SKY_MAP_MEMORY.flow,...
                                                  SKY_MAP_MEMORY.deltaF,...
                                                  SKY_MAP_MEMORY.numFreqs,...
                                                  GPStoGreenwichMeanSiderealTime(GPSstart),...
                                                  SKY_MAP_MEMORY.glm,...
                                                  SKY_MAP_MEMORY.g1lm,...
                                                  SKY_MAP_MEMORY.g2lm,...
                                                  SKY_MAP_MEMORY.det1,...
                                                  SKY_MAP_MEMORY.det2,...
                                                  SKY_MAP_MEMORY.isSpH,...
                                                  SKY_MAP_MEMORY.coord1,...
                                                  SKY_MAP_MEMORY.coord2,...
                                                  SKY_MAP_MEMORY.map);
  P1 = P1 + simP1;
  P2 = P2 + simP2;
  C = C + simC;
end
           
vSpH = doSpH(vSpH,params,C,P1,P2,GPSstart);

