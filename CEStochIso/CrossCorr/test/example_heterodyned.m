% Script which tests stochastic pipeline using heterodyned data.
% It uses the same data as example_params.txt but first converts
% the frame data to heterodyned matlab files.
% The results for the theoretical sigma are virtually identical
% to that obtained using example_params.txt, as expected.
%
clear all;

stochastic_paths;

fbase = 100;
phi = 0;

% Flag for matlab output (frame output does not seem to
% work, it produces zeroes for the data)
matlabOut = true;

paramsFile = 'example_params.txt';
jobsFile = 'example_jobs.txt';
jobNumber = 2;

params = setStochasticParams(paramsFile, jobsFile, jobNumber);
%checkParamsStochastic(params);
params = loadAuxiliaryInput(params);

tmpDuration = 6*params.segmentDuration;
dataStartTime = 1126627298 + 2 - params.bufferSecs1;

fprintf('Heterodyning H1 data...\n');
[adcdata, dataOK] = readTimeSeriesData(params.channelName1, ...
                                       dataStartTime, ...
                                       tmpDuration, ...
                                       params.frameDuration1,...
                                       params.gpsTimesFile1, ...
                                       params.frameCacheFile1);
q = round(1/(adcdata.deltaT*params.resampleRate1));
adcdata.data = resample(adcdata.data, 1, q, params.nResample1, params.betaParam1);
adcdata.deltaT = q*adcdata.deltaT;
if (params.doHighPass1)
  adcdata.data = cascadefilter(adcdata.data, params.highPassOrder1, params.highPassFreq1, params.resampleRate1);
end;
adcdata.fbase = fbase;
adcdata.phase = phi;
t = adcdata.deltaT*[0:1:length(adcdata.data)-1].';
adcdata.data = adcdata.data.*exp(-i*2*pi*adcdata.fbase.*t + adcdata.phase);

if (matlabOut)
  fname = [ params.ifo1(1) '-' params.ifo1 '_HETERODYNED' '-' num2str(dataStartTime) '-' num2str(tmpDuration) '.mat' ];
  save([ 'old/output/' fname ], 'adcdata');
else
  data.channel = params.channelName1;
  data.data = adcdata.data;
  data.type = 'dc';
  data.mode = 'a';
  fname = [ params.ifo1(1) '-' params.ifo1 '_HETERODYNED' '-' num2str(dataStartTime) '-' num2str(tmpDuration) '.gwf' ];
  mkframe([ 'old/output/' fname ], data, 'n', tmpDuration, dataStartTime);
end;

f = fopen('old/output/frameFilesH.2.txt', 'w');
fprintf(f, [ pwd '/old/output/' fname '\n' ]);
fclose(f);

fprintf('Heterodyning L1 data...\n');
[adcdata, dataOK] = readTimeSeriesData(params.channelName2, ...
                                       dataStartTime, ...
                                       tmpDuration, ...
                                       params.frameDuration2,...
                                       params.gpsTimesFile2, ...
                                       params.frameCacheFile2);

q = round(1/(adcdata.deltaT*params.resampleRate2));
adcdata.data = resample(adcdata.data, 1, q, params.nResample2, params.betaParam2);
adcdata.deltaT = q*adcdata.deltaT;
if (params.doHighPass2)
  adcdata.data = cascadefilter(adcdata.data, params.highPassOrder2, params.highPassFreq2, params.resampleRate2);
end;
adcdata.fbase = fbase;
adcdata.phase = phi;
t = adcdata.deltaT*[0:1:length(adcdata.data)-1].';
adcdata.data = adcdata.data.*exp(-i*2*pi*adcdata.fbase.*t + adcdata.phase);

if (matlabOut)
  fname = [ params.ifo2(1) '-' params.ifo2 '_HETERODYNED' '-' num2str(dataStartTime) '-' num2str(tmpDuration) '.mat' ];
  save([ 'old/output/' fname ], 'adcdata');
else
  % Frame file
  data.channel = params.channelName2;
  data.data = adcdata.data;
  data.type = 'dc';
  data.mode = 'a';

  fname = [ params.ifo2(1) '-' params.ifo2 '_HETERODYNED' '-' num2str(dataStartTime) '-' num2str(tmpDuration) '.gwf' ];
  mkframe([ 'old/output/' fname ], data, 'n', tmpDuration, dataStartTime);
end;

f = fopen('old/output/frameFilesL.2.txt', 'w');
fprintf(f, [ pwd '/old/output/' fname '\n' ]);
fclose(f);

delete('example_heterodyned.out');
diary('example_heterodyned.out');
stochastic(paramsFile, jobsFile, jobNumber);
diary off;
