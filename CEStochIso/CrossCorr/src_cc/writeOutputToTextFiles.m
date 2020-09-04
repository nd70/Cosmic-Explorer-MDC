function writeOutputToTextFiles(params, I, K, segmentStartTime, result)
%
% writes the data to text files
%
% $Id: writeToOutputFiles.m,v 1.2 2007-07-18 01:18:34 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if params.writeStatsToFiles
    if params.doDirectional
      % write whole ccStat time series to fiule
      numTimes=length(result.ccStat.data);
      fprintf(params.ccstats_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
              [segmentStartTime*transpose(ones(numTimes,1)); ...
               sqrt(result.ccVar)*transpose(ones(numTimes,1)); ...
               (-numTimes/2:numTimes/2-1)*result.ccStat.deltaT; ...
               transpose((result.ccStat.data))]);
    else
      % write cc stats, etc. to output files 
      if params.heterodyned
        fprintf(params.ccstats_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
                segmentStartTime, real(result.ccStat), imag(result.ccStat), sqrt(result.ccVar));
      else 
        fprintf(params.ccstats_fid(K), '%9.1f\t%e\t%e\n', ...
                segmentStartTime, result.ccStat, sqrt(result.ccVar));
      end
    end
    if params.doAllSkyComparison
      fprintf(params.ccstatsAllSky_fid(K), '%9.1f\t%e\t%e\n', ...
              segmentStartTime, result.ccStatAllSky, sqrt(result.ccVarAllSky));
    end
  end

  if params.writeNaiveSigmasToFiles
    % write "naive" sigmas, etc. to output files 
    fprintf(params.naivesigmas_fid(K), '%9.1f\t%e\t%e\n', ...
            segmentStartTime, sqrt(result.naiVar), sqrt(result.ccVar));
    if params.doAllSkyComparison
      fprintf(params.naivesigmasAllSky_fid(K), '%9.1f\t%e\t%e\n', ...
              segmentStartTime, sqrt(result.naiVarAllSky), sqrt(result.ccVarAllSky));
    end
  end

  if params.writeSpectraToFiles
    % write cc spectra to output files
    fprintf(params.ccspectra_fid(K), '%9.1f\t%e\t%e\t%e\t%e\n', ...
            [segmentStartTime*transpose(ones(params.numFreqs,1)); ...
             sqrt(result.ccVar)*transpose(ones(params.numFreqs,1)); transpose(params.f); ...
             transpose(real(result.ccSpec.data)); transpose(imag(result.ccSpec.data))]);
    if params.doAllSkyComparison
      fprintf(params.ccspectraAllSky_fid(K), '%9.1f\t%e\t%e\t%e\t%e\n', ...
              [segmentStartTime*transpose(ones(params.numFreqs,1)); ...
               sqrt(result.ccVarAllSky)*transpose(ones(params.numFreqs,1)); transpose(params.f); ...
               transpose(real(result.ccSpecAllSky.data)); transpose(imag(result.ccSpecAllSky.data))]);
    end
  end

  if params.writeSensIntsToFiles
    % write sensitivity integrands to output files
    fprintf(params.sensints_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
            [segmentStartTime*ones(params.numFreqs,1)'; ...
             sqrt(result.ccVar)*transpose(ones(params.numFreqs,1)); ...
             transpose(params.f); transpose(result.sensInt.data)]);
    if params.doAllSkyComparison
      fprintf(params.sensintsAllSky_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
              [segmentStartTime*ones(params.numFreqs,1)'; ...
               sqrt(result.ccVarAllSky)*transpose(ones(params.numFreqs,1)); ...
               transpose(params.f); transpose(result.sensIntAllSky.data)]);
    end
  end

  
  if params.writeOptimalFiltersToFiles
    % write optimal filters to output files
    fprintf(params.optimal_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
            [segmentStartTime*ones(params.numFreqs,1)'; ...
             transpose(params.f); transpose(real(result.Q.data));transpose(imag(result.Q.data))]);
    if params.doAllSkyComparison
      fprintf(params.optimalAllSky_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
              [segmentStartTime*ones(params.numFreqs,1)'; ...
               transpose(params.f); transpose(real(result.QAllSky.data));transpose(imag(result.QAllSky.data))]);
    end
  end

  if params.writeOverlapReductionFunctionToFiles
    % write overlap reduction function to output files
    fprintf(params.orf_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
            [segmentStartTime*ones(params.numFreqs,1)'; ...
             transpose(params.f); transpose(real(params.gamma.data));transpose(imag(params.gamma.data))]);
  end

  if params.writeCalPSD1sToFiles
    % write calibrated power spectra to output files
    fprintf(params.calpsd1_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
            [segmentStartTime*ones(params.numFreqs,1)'; ...
             sqrt(result.ccVar)*ones(params.numFreqs,1)'; ...
             transpose(params.f); transpose(result.calPSD1_avg.data)]);
  end

  if params.writeCalPSD2sToFiles
    % write calibrated power spectra to output files
    fprintf(params.calpsd2_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
            [segmentStartTime*ones(params.numFreqs,1)'; ...
             sqrt(result.ccVar)*ones(params.numFreqs,1)'; ...
             transpose(params.f); transpose(result.calPSD2_avg.data)]);
  end
  
  % Write out bad GPS times here unless combining data. If combined, the bad GPS times will be written out
  % by finishCombine
  if (params.doBadGPSTimes & ~params.doCombine)
    if (sqrt(result.ccVar/result.naiVar) > params.maxDSigRatio | sqrt(result.ccVar/result.naiVar) < params.minDSigRatio)
      fprintf(params.badGPSTimes_fid(K),'%9.1f\n', segmentStartTime);
    end;
  end;

return;
