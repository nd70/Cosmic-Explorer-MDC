function writeOutputToMatFile(params, I, K, segmentStartTime, result)
%
% Writes the results for segment I trial K to the matfile.
%
% This routine uses the matfile object to apend the results from segment I to
% an existing Matlab file. The matfile object is able to read and write
% partial data in Matlab save files without reading the entire file, so it
% is very efficient.
%
% The resulting matfile contains the outputs from each segment I and trial K
% arranged in a matrix result(I, K).
%
% $Id: writeToOutputFiles.m,v 1.2 2007-07-18 01:18:34 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % When writing to matfile we will store the results as an array of structs
  % containing the output, indexed by segment. Since we only want to save select
  % data we will create a local result and copy the desired output into it
  % before saving.
  % No trial dimension is used because segmentStartTime must be the same for
  % all trials.
  params.matFile.segmentStartTime(I, 1) = segmentStartTime;

  % ccStat
  if params.writeStatsToFiles
    params.matFile.ccStat(I, K) = result.ccStat;
  end

  % ccVar (written in many cases)
  % When writing text files, standard deviations are written rather than
  % variances, and the post-processing scripts expect to see standard
  % deviations, so we follow the same convention for matfile output.
  if (params.writeStatsToFiles ...
      | params.writeNaiveSigmasToFiles ...
      | params.writeSpectraToFiles ...
      | params.writeSensIntsToFiles ...
      | params.writeCalPSD1sToFiles ...
      | params.writeCalPSD2sToFiles)
    params.matFile.ccSigma(I, K) = sqrt(result.ccVar);
  end

  % ccStatAllSky
  if (params.writeStatsToFiles & params.doAllSkyComparison)
    params.matFile.ccStatAllSky(I, K) = result.ccStatAllSky;
  end

  % ccVarAllSky (written in many cases)
  if (params.doAllSkyComparison ...
      & (params.writeStatsToFiles ...
       | params.writeNaiveSigmasToFiles ...
       | params.writeSpectraToFiles ...
       | params.writeSensIntsToFiles ...
       | params.writeCalPSD1sToFiles ...
       | params.writeCalPSD2sToFiles))
    params.matFile.ccSigmaAllSky(I, K) = sqrt(result.ccVarAllSky);
  end

  % naiVar
  if params.writeNaiveSigmasToFiles
    % write "naive" sigmas, etc. to output files 
    params.matFile.naiSigma(I, K) = sqrt(result.naiVar);
  end

  % naiVarAllSky
  if (params.doAllSkyComparison & params.writeNaiveSigmasToFiles)
    params.matFile.naiSigmaAllSky(I, K) = sqrt(result.naiVarAllSky);
  end

  % ccSpec
  if params.writeSpectraToFiles
    params.matFile.ccSpec(I, K) = result.ccSpec;
  end

  % ccSpecAllSky
  if (params.writeSpectraToFiles & params.doAllSkyComparison)
    params.matFile.ccSpecAllSky(I, K) = result.ccSpecAllSky;
  end

  % sensInt
  if params.writeSensIntsToFiles
    % write sensitivity integrands to output files
    params.matFile.sensInt(I, K) = result.sensInt;
  end

  % sensIntAllSky
  if (params.doAllSkyComparison & params.writeSensIntsToFiles)
    % write sensitivity integrands to output files
    params.matFile.sensIntAllSky(I, K) = result.sensIntAllSky;
  end

  % Q
  if params.writeOptimalFiltersToFiles
    % write optimal filters to output files
    params.matFile.Q(I, K) = result.Q;
  end

  % QAllSky
  if (params.doAllSkyComparison & params.writeOptimalFiltersToFiles)
    % write optimal filters to output files
    params.matFile.QAllSky(I, K) = result.QAllSky;
  end

  % calPSD1_avg
  if params.writeCalPSD1sToFiles
    % write calibrated power spectra to output files
    params.matFile.calPSD1_avg(I, K) = result.calPSD1_avg;
  end
  
  % calPSD2_avg
  if params.writeCalPSD2sToFiles
    % write calibrated power spectra to output files
    params.matFile.calPSD2_avg(I, K) = result.calPSD2_avg;
  end

  % badGPSTimes
  %
  % Unlike the other data products, badGPSTimes is only updated when the segment is bad, not for every segment.
  % It is stored as a cellarray of vectors, where the cellarray is indexed by trial K ie. badGPSTimes{K} is the
  % vector of bad GPS times in trial K. This is also problematic for writing using matfile so we need to
  % instead load the list of bad times, append the new one it and save the new version, however this is infrequent.
  if params.doBadGPSTimes
    if (sqrt(result.ccVar/result.naiVar) > params.maxDSigRatio | sqrt(result.ccVar/result.naiVar) < params.minDSigRatio)
      badGPSTimes = params.matFile.badGPSTimes;
      badGPSTimes{K}(end+1, 1) = segmentStartTime;
      params.matFile.badGPSTimes = badGPSTimes;
    end;
  end;

return;
