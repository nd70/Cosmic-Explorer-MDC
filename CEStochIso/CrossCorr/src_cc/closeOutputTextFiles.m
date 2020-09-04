function params=closeOutputTextFiles(params)

% closes all output text files

% input and output: params struct
%
% $Id: closeOutputFiles.m,v 1.1 2007-06-26 00:18:48 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % CLOSE FILES USED TO SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if params.doDirectional
    try
      Sky = params.Sky;
    catch
      Sky = [];
    end;
    if length(Sky)==0
      Sky=[]; % this is for backward compatibility reasons
    end
    try
      save(params.ccStatSkySetCurrentFilename,'Sky');
    catch
      error('Error saving ccStatSkySet file');
    end;
    clear Sky;
    params.Sky={};
    try
      Sky = params.metaSky;
    catch
      Sky = [];
    end;
    if length(Sky)==0
      Sky=[]; % this is for backward compatibility reasons
    end
    try
      save(params.ccStatSkyfilename, 'Sky');
    catch
      error('Error saving ccStatSky file');
    end;
    clear Sky;
    params.metaSky={};
  end;

  for K=1:params.numTrials

    if params.doCombine
      fclose(params.combined_ccstats_fid(K));
      fclose(params.combined_ccspectra_fid(K));
      fclose(params.combined_sensints_fid(K));
      if params.doBadGPSTimes
        fclose(params.badGPSTimes_fid(K));
      end
    end

    if params.writeStatsToFiles
      fclose(params.ccstats_fid(K));
      if params.doAllSkyComparison
        fclose(params.ccstatsAllSky_fid(K));
      end
    end

    if params.writeNaiveSigmasToFiles
      fclose(params.naivesigmas_fid(K));
      if params.doAllSkyComparison
        fclose(params.naivesigmasAllSky_fid(K));
      end
    end

    if params.writeSpectraToFiles
      fclose(params.ccspectra_fid(K));
      if params.doAllSkyComparison
        fclose(params.ccspectraAllSky_fid(K));
      end
    end

    if params.writeSensIntsToFiles
      fclose(params.sensints_fid(K));
      if params.doAllSkyComparison
        fclose(params.sensintsAllSky_fid(K));
      end
    end

    if params.writeCoherenceToFiles
      fclose(params.coherence_fid(K));
    end
    
    if params.writeOptimalFiltersToFiles
      fclose(params.optimal_fid(K));
      if params.doAllSkyComparison
        fclose(params.optimalAllSky_fid(K));
      end
    end

    if params.writeOverlapReductionFunctionToFiles
      fclose(params.orf_fid(K));
    end

    if params.writeCalPSD1sToFiles
      fclose(params.calpsd1_fid(K));
    end

    if params.writeCalPSD2sToFiles
      fclose(params.calpsd2_fid(K));
    end

  end

return;