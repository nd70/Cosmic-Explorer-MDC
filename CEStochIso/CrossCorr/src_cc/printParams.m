function printParams(params)

% print a copy of the params struct into the output directory for
% future reference.
% 
% input: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer
% sballmer@caltech.edu 
%
% $Id: printParams.m,v 1.8 2009-02-27 19:56:54 ethrane Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% write date/time, params, etc. to a file
filename = [params.outputFilePrefix '_params.job' num2str(params.jobNumber) '.dat'];
fid = fopen(filename, 'w');
if (fid == -1)
  error(sprintf('Cannot open file \''%s\'' for writing.', filename));
end;

fprintf(fid, '%% Date and time of this run: %s\n', params.ddmmyyyyhhmmss);
if params.UnphysicalTimeShift~=0
  fprintf(fid, '%% Warning: params.UnphysicalTimeShift is not 0, end result will be garbage!\n');
end
fprintf(fid, '%s\n', ...
        '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid, '%% matapps CVS tag name string = %s\n', '$Name: not supported by cvs2svn $');
fprintf(fid, '%% stochastic.m CVS version string = %s\n', '$Id: printParams.m,v 1.8 2009-02-27 19:56:54 ethrane Exp $');
fprintf(fid, '%% paramsFile = %s\n', params.paramsFile);
fprintf(fid, '%% jobsFile = %s\n', params.jobsFile);
fprintf(fid, '%% jobNumber = %d\n', params.jobNumber);
fprintf(fid, '%% doDirectional = %d\n', params.doDirectional);
fprintf(fid, '%% doSphericalHarmonics = %d\n', params.doSphericalHarmonics);
fprintf(fid, '%% doNarrowbandRadiometer = %d\n', params.doNarrowbandRadiometer);
fprintf(fid, '%% doAllSkyComparison = %d\n', params.doAllSkyComparison);
fprintf(fid, '%% doFreqMask = %d\n', params.doFreqMask);
fprintf(fid, '%% doHighPass1 = %d\n', params.doHighPass1);
fprintf(fid, '%% doHighPass2 = %d\n', params.doHighPass1);
fprintf(fid, '%% doMonteCarlo = %d\n', params.doMonteCarlo);
fprintf(fid, '%% doSimulatedPointSource = %d\n', params.doSimulatedPointSource);
fprintf(fid, '%% doSimulatedSkyMap = %d\n', params.doSimulatedSkyMap);
fprintf(fid, '%% doMCoffset = %d\n', params.doMCoffset);
fprintf(fid, '%% doConstTimeShift = %d\n', params.doConstTimeShift);
fprintf(fid, '%% doOverlap = %d\n', params.doOverlap);
fprintf(fid, '%% doCombine = %d\n', params.doCombine);
fprintf(fid, '%% doBadGPSTimes = %d\n', params.doBadGPSTimes);
fprintf(fid, '%% badGPSTimesFile = %s\n', params.badGPSTimesFile); %Feb27 cet
fprintf(fid, '%% doShift1 = %d\n', params.doShift1);
fprintf(fid, '%% doShift2 = %d\n', params.doShift2);
fprintf(fid, '%% doInjFromFile1 = %d\n', params.doInjFromFile1);
fprintf(fid, '%% doInjFromFile2 = %d\n', params.doInjFromFile2);
fprintf(fid, '%% heterodyned = %d\n', params.heterodyned);
fprintf(fid, '%% useDatafindServer = %d\n', params.useDatafindServer);
fprintf(fid, '%% writeResultsToScreen = %d\n', params.writeResultsToScreen);
fprintf(fid, '%% writeStatsToFiles = %d\n',params.writeStatsToFiles);
fprintf(fid, '%% writeOutputToMatFile = %d\n',params.writeOutputToMatFile);
fprintf(fid, '%% writeNaiveSigmasToFiles = %d\n',params.writeNaiveSigmasToFiles);
fprintf(fid, '%% writeSpectraToFiles = %d\n', params.writeSpectraToFiles);
fprintf(fid, '%% writeSensIntsToFiles = %d\n', params.writeSensIntsToFiles);
fprintf(fid, '%% writeCoherenceToFiles = %d\n', params.writeCoherenceToFiles);
fprintf(fid, '%% writeCohFToFiles = %d\n', params.writeCohFToFiles);
fprintf(fid, '%% writeOptimalFiltersToFiles = %d\n',params.writeOptimalFiltersToFiles);
fprintf(fid, '%% writeOverlapReductionFunctionToFiles = %d\n',params.writeOverlapReductionFunctionToFiles);
fprintf(fid, '%% writeCalPSD1sToFiles = %d\n', params.writeCalPSD1sToFiles);
fprintf(fid, '%% writeCalPSD2sToFiles = %d\n', params.writeCalPSD2sToFiles);
fprintf(fid, '%% ifo1 = %s\n', params.ifo1);
fprintf(fid, '%% ifo2 = %s\n', params.ifo2);
if (isnan(params.azimuth1))
  fprintf(fid, '%% azimuth1 NOT SET\n', params.azimuth1);
else
  fprintf(fid, '%% azimuth1 = %f degrees\n', params.azimuth1);
end;
if (isnan(params.azimuth2))
  fprintf(fid, '%% azimuth2 NOT SET\n', params.azimuth2);
else
  fprintf(fid, '%% azimuth2 = %f degrees\n', params.azimuth2);
end;
fprintf(fid, '%% segmentDuration = %f sec\n', params.segmentDuration);
fprintf(fid, '%% numSegmentsPerInterval = %d\n', params.numSegmentsPerInterval);
fprintf(fid, '%% ignoreMidSegment = %d\n', params.ignoreMidSegment);
fprintf(fid, '%% deltaF = %f Hz\n', params.deltaF);
fprintf(fid, '%% flow = %f Hz\n', params.flow);
fprintf(fid, '%% fhigh = %f Hz\n', params.fhigh);
if length(params.alphaExp)==0
  fprintf(fid, '%% alphaExp = []\n');
else
  fprintf(fid, '%% alphaExp = %f\n', params.alphaExp);
end
if isstruct(params.fRef)
  fprintf(fid, '%% fRef.data size is %d \n%% fRef.flow = %f Hz\n%% fRef.deltaF = %f Hz\n%% fRef.symmetry = %d\n',...
  length(params.fRef.data), params.fRef.flow,params.fRef.deltaF,params.fRef.symmetry);
else
  fprintf(fid, '%% fRef = %f Hz\n', params.fRef);
end
fprintf(fid, '%% maxSegmentsPerMatfile = %f Hz\n', params.maxSegmentsPerMatfile);
fprintf(fid, '%% useSignalSpectrumHfFromFile = %d\n', params.useSignalSpectrumHfFromFile);
fprintf(fid, '%% HfFile = %s\n', params.HfFile);
fprintf(fid, '%% HfFileInterpolateLogarithmic = %f\n', params.HfFileInterpolateLogarithmic);
fprintf(fid, '%% SpHLmax = %d\n', params.SpHLmax);
fprintf(fid, '%% useSkyPatternFile = %f\n', params.useSkyPatternFile);
fprintf(fid, '%% SkyPatternFile = %s\n', params.SkyPatternFile);
fprintf(fid, '%% SkyPatternRightAscensionNumPoints = %f\n', params.SkyPatternRightAscensionNumPoints);
fprintf(fid, '%% SkyPatternDeclinationNumPoints = %f\n', params.SkyPatternDeclinationNumPoints);
fprintf(fid, '%% maxCorrelationTimeShift = %f\n', params.maxCorrelationTimeShift);
fprintf(fid, '%% UnphysicalTimeShift = %f\n', params.UnphysicalTimeShift);
fprintf(fid, '%% resampleRate1 = %f Hz\n', params.resampleRate1);
fprintf(fid, '%% resampleRate2 = %f Hz\n', params.resampleRate2);
if (isnan(params.fbase1))
  fprintf(fid, '%% fbase1 NOT SET\n', params.fbase1);
else
  fprintf(fid, '%% fbase1 = %f Hz\n', params.fbase1);
end;
if (isnan(params.fbase2))
  fprintf(fid, '%% fbase2 NOT SET\n', params.fbase2);
else
  fprintf(fid, '%% fbase2 = %f Hz\n', params.fbase2);
end;
fprintf(fid, '%% bufferSecs1 = %f sec\n', params.bufferSecs1);
fprintf(fid, '%% bufferSecs2 = %f sec\n', params.bufferSecs2);
fprintf(fid, '%% ASQchannel1 = %s\n', params.ASQchannel1);
fprintf(fid, '%% ASQchannel2 = %s\n', params.ASQchannel2);
fprintf(fid, '%% injChannel1 = %s\n', params.injChannel1);
fprintf(fid, '%% injChannel2 = %s\n', params.injChannel2);
fprintf(fid, '%% injPrefix1 = %s\n', params.injPrefix1);
fprintf(fid, '%% injPrefix2 = %s\n', params.injPrefix2);
fprintf(fid, '%% injScale1 = %f\n', params.injScale1);
fprintf(fid, '%% injScale2 = %f\n', params.injScale2);
fprintf(fid, '%% frameType1 = %s\n', params.frameType1);
fprintf(fid, '%% frameType2 = %s\n', params.frameType2);
fprintf(fid, '%% frameDuration1 = %f sec\n', params.frameDuration1);
fprintf(fid, '%% frameDuration2 = %f sec\n', params.frameDuration2);
fprintf(fid, '%% injFrameDuration1 = %f sec\n', params.injFrameDuration1);
fprintf(fid, '%% injFrameDuration2 = %f sec\n', params.injFrameDuration2);
fprintf(fid, '%% hannDuration1 = %f sec\n', params.hannDuration1);
fprintf(fid, '%% hannDuration2 = %f sec\n', params.hannDuration2);
fprintf(fid, '%% nResample1 = %d \n', params.nResample1);
fprintf(fid, '%% nResample2 = %d \n', params.nResample2);
fprintf(fid, '%% betaParam1 = %d\n', params.betaParam1);
fprintf(fid, '%% betaParam2 = %d\n', params.betaParam2);
fprintf(fid, '%% highPassFreq1 = %f Hz\n', params.highPassFreq1);
fprintf(fid, '%% highPassFreq2 = %f Hz\n', params.highPassFreq2);
fprintf(fid, '%% highPassOrder1 = %d\n', params.highPassOrder1);
fprintf(fid, '%% highPassOrder2 = %d\n', params.highPassOrder2);
fprintf(fid, '%% freqsToRemove = %f\n', params.freqsToRemove);
fprintf(fid, '%% nBinsToRemove = %d\n', params.nBinsToRemove);
fprintf(fid, '%% numTrials = %d\n', params.numTrials);
fprintf(fid, '%% signalType = %s\n', params.signalType);
fprintf(fid, '%% simOmegaRef1 = %f\n', params.simOmegaRef1);
fprintf(fid, '%% simOmegaRef2 = %f\n', params.simOmegaRef2);
fprintf(fid, '%% simulationPath = %s\n', params.simulationPath);
fprintf(fid, '%% simulatedPointSourcesFile = %s\n', params.simulatedPointSourcesFile);
fprintf(fid, '%% simulatedPointSourcesPowerSpec = %s\n', params.simulatedPointSourcesPowerSpec);
fprintf(fid, '%% simulatedPointSourcesInterpolateLogarithmic = %d\n', params.simulatedPointSourcesInterpolateLogarithmic);
fprintf(fid, '%% simulatedPointSourcesBufferDepth = %f sec\n', params.simulatedPointSourcesBufferDepth);
fprintf(fid, '%% simulatedPointSourcesHalfRefillLength = %f sec\n', params.simulatedPointSourcesHalfRefillLength);
fprintf(fid, '%% simulatedPointSourcesNoRealData = %d\n', params.simulatedPointSourcesNoRealData);
fprintf(fid, '%% simulatedPointSourcesMakeIncoherent = %d\n', params.simulatedPointSourcesMakeIncoherent);
fprintf(fid, '%% simulatedSkyMapFile = %s\n', params.simulatedSkyMapFile);
fprintf(fid, '%% simulatedSkyMapFileType = %d\n', params.simulatedSkyMapFileType);
fprintf(fid, '%% simulatedSkyMapFileNumBins = %d\n', params.simulatedSkyMapFileNumBins);
fprintf(fid, '%% simulatedSkyMapInjectAsSpH = %d\n', params.simulatedSkyMapInjectAsSpH);
fprintf(fid, '%% simulatedSkyMapConvertLmax = %d\n', params.simulatedSkyMapConvertLmax);
fprintf(fid, '%% simulatedSkyMapConvertDeg = %d\n', params.simulatedSkyMapConvertDeg);
fprintf(fid, '%% simulatedSkyMapInjectTimeDomain = %d\n', params.simulatedSkyMapInjectTimeDomain);
fprintf(fid, '%% simulatedSkyMapPowerSpec = %s\n', params.simulatedSkyMapPowerSpec);
fprintf(fid, '%% simulatedSkyMapInterpolateLogarithmic = %d\n', params.simulatedSkyMapInterpolateLogarithmic);
fprintf(fid, '%% simulatedSkyMapBufferDepth = %f sec\n', params.simulatedSkyMapBufferDepth);
fprintf(fid, '%% simulatedSkyMapHalfRefillLength = %f sec\n', params.simulatedSkyMapHalfRefillLength);
fprintf(fid, '%% simulatedSkyMapNoRealData = %d\n', params.simulatedSkyMapNoRealData);
fprintf(fid, '%% simulatedSkyMapMakeIncoherent = %d\n', params.simulatedSkyMapMakeIncoherent);
fprintf(fid, '%% alphaBetaFile1 = %s\n', params.alphaBetaFile1);
fprintf(fid, '%% minMCoff = %f\n', params.minMCoff);
fprintf(fid, '%% maxMCoff = %f\n', params.maxMCoff);
fprintf(fid, '%% ConstTimeShift = %f\n', params.ConstTimeShift);
fprintf(fid, '%% ShiftTime1 = %f\n', params.ShiftTime1);
fprintf(fid, '%% ShiftTime2 = %f\n', params.ShiftTime2);
fprintf(fid, '%% maxDSigRatio = %f\n', params.maxDSigRatio);
fprintf(fid, '%% minDSigRatio = %f\n', params.minDSigRatio);
fprintf(fid, '%% minDataLoadLength = %f\n', params.minDataLoadLength);
fprintf(fid, '%% alphaBetaFile1 = %s\n', params.alphaBetaFile1);
fprintf(fid, '%% alphaBetaFile2 = %s\n', params.alphaBetaFile2);
fprintf(fid, '%% calCavGainFile1 = %s\n', params.calCavGainFile1);
fprintf(fid, '%% calCavGainFile2 = %s\n', params.calCavGainFile2);
fprintf(fid, '%% calResponseFile1 = %s\n', params.calResponseFile1);
fprintf(fid, '%% calResponseFile2 = %s\n', params.calResponseFile2);
if ~params.intermediate %cet--------------------------------------------------
  fprintf(fid, '%% gpsTimesPath1 = %s\n', params.gpsTimesPath1);
  fprintf(fid, '%% gpsTimesPath2 = %s\n', params.gpsTimesPath2);
  fprintf(fid, '%% frameCachePath1 = %s\n', params.frameCachePath1);
  fprintf(fid, '%% frameCachePath2 = %s\n', params.frameCachePath2);
else
  fprintf(fid, '%% intFrameCachePath = %s\n', params.intFrameCachePath);
end %-------------------------------------------------------------------------
fprintf(fid, '%% injGPSTimesPath1 = %s\n', params.injGPSTimesPath1);
fprintf(fid, '%% injGPSTimesPath2 = %s\n', params.injGPSTimesPath2);
fprintf(fid, '%% injFrameCachePath1 = %s\n', params.injFrameCachePath1);
fprintf(fid, '%% injFrameCachePath2 = %s\n', params.injFrameCachePath2);
fprintf(fid, '%% outputFilePrefix = %s\n', params.outputFilePrefix);
  
fclose(fid);
