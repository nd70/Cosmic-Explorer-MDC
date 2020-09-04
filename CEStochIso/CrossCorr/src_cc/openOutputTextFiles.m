function params=openOutputTextFiles(params)

% opens all output files for writing when using text format
% 
% input and output: params struct
%
% $Id: openOutputFiles.m,v 1.1 2007-06-26 00:18:48 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % OPEN FILES FOR SAVING RESULTS (if desired) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if params.doDirectional
    % open files for cc stats, etc. and write header info
    params.ccStatSkyfilename  = [params.outputFilePrefix '_ccstatsSky.job' ...
                        num2str(params.jobNumber) '.trial' num2str(1) '.mat'];
    params.ccStatSkySetPrefix = [params.outputFilePrefix '_ccstatsSkySet'];
    params.ccStatSkySetSuffix = ['.job' num2str(params.jobNumber) '.trial' num2str(1) '.mat'];
    params.skyIndex=1;
    params.metaSky={};
    params.skySetNumber=1;
    params.skySetSegmentOffset=0;
    params.ccStatSkySetCurrentFilename= [params.ccStatSkySetPrefix num2str(params.skySetNumber) params.ccStatSkySetSuffix];
    params.Sky={};
  end

  for K=1:params.numTrials

    if params.doCombine
      % open files for combined results
      filename = [params.outputFilePrefix '_combined_ccstats.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.combined_ccstats_fid(K) = fopen(filename, 'w');
      fprintf(params.combined_ccstats_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.combined_ccstats_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.combined_ccstats_fid(K), '%s\t%s\t%s\n', ...
              '%start sec', 'CC statistic', 'theor sigma');

      filename = [params.outputFilePrefix '_combined_ccspectra.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.combined_ccspectra_fid(K) = fopen(filename, 'w');
      fprintf(params.combined_ccspectra_fid(K), '%%Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.combined_ccspectra_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.combined_ccspectra_fid(K), '%s\t%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', ...
              'cc spec (real)', 'cc spec (imag)');

      filename = [params.outputFilePrefix '_combined_sensints.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.combined_sensints_fid(K) = fopen(filename, 'w');
      fprintf(params.combined_sensints_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.combined_sensints_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.combined_sensints_fid(K), '%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', ...
              'sens integrand (1/sec)');

    end

    if params.writeStatsToFiles
      % open files for cc stats, etc. and write header info
      filename = [params.outputFilePrefix '_ccstats.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.ccstats_fid(K) = fopen(filename, 'w');

      fprintf(params.ccstats_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.ccstats_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      if params.doDirectional
        fprintf(params.ccstats_fid(K), '%s\t%s\t%s\t%s\n', ...
                '%start sec', 'theor sigma', 'shift(sec)', ...
                'CC stat');
      else
        if params.heterodyned
          fprintf(params.ccstats_fid(K), '%s\t%s\t%s\t%s\n', ...
                  '%start sec', 'Re(CC stat)', 'Im(CC stat)', ...
                  'theor sigma');
        else
          fprintf(params.ccstats_fid(K), '%s\t%s\t%s\n', ...
                  '%start sec', 'CC statistic', 'theor sigma');
        end
      end
      if params.doAllSkyComparison
        % open files for cc stats, etc. and write header info
        filename = [params.outputFilePrefix '_AllSky_ccstats.job' ...
                    num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
        params.ccstatsAllSky_fid(K) = fopen(filename, 'w');

        fprintf(params.ccstatsAllSky_fid(K), '%% Date and time of this run: %s\n', ...
                params.ddmmyyyyhhmmss);
        fprintf(params.ccstatsAllSky_fid(K), '%s\n', ...
                '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(params.ccstatsAllSky_fid(K), '%s\t%s\t%s\n', ...
                '%start sec', 'CC statistic', 'theor sigma');
      end
    end

    % This has been moved from within the doCombine test so that badGPSTimes can be written
    % for both combined and uncombined data
    if params.doBadGPSTimes
      filename = [params.outputFilePrefix '_badGPSTimes.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.badGPSTimes_fid(K) = fopen(filename, 'w');
      fprintf(params.badGPSTimes_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.badGPSTimes_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    end;
    
    if params.writeNaiveSigmasToFiles
      % open files for naive sigmas (those calculated w/o the sliding
      % PSD estimator) , etc. and write header info
      filename = [params.outputFilePrefix '_naivesigmas.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.naivesigmas_fid(K) = fopen(filename, 'w');

      fprintf(params.naivesigmas_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.naivesigmas_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.naivesigmas_fid(K), '%s\t%s\t%s\n', ...
              '%start sec', 'naive theor sigma', 'theor sigma');
      if params.doAllSkyComparison
        % open files for naive sigmas (those calculated w/o the sliding
        % PSD estimator) , etc. and write header info
        filename = [params.outputFilePrefix '_AllSky_naivesigmas.job' ...
                    num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
        params.naivesigmasAllSky_fid(K) = fopen(filename, 'w');

        fprintf(params.naivesigmasAllSky_fid(K), '%% Date and time of this run: %s\n', ...
                params.ddmmyyyyhhmmss);
        fprintf(params.naivesigmasAllSky_fid(K), '%s\n', ...
                '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(params.naivesigmasAllSky_fid(K), '%s\t%s\t%s\n', ...
                '%start sec', 'naive theor sigma', 'theor sigma');
      end
    end;

    if params.writeSpectraToFiles
      % open files for cc spectra, etc. and write header info
      filename = [params.outputFilePrefix '_ccspectra.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.ccspectra_fid(K) = fopen(filename, 'w');

      fprintf(params.ccspectra_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.ccspectra_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.ccspectra_fid(K), '%s\t%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', ...
              'cc spec (real)', 'cc spec (imag)');

      if params.doAllSkyComparison
        % open files for cc spectra, etc. and write header info
        filename = [params.outputFilePrefix '_AllSky_ccspectra.job' ...
                    num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
        params.ccspectraAllSky_fid(K) = fopen(filename, 'w');

        fprintf(params.ccspectraAllSky_fid(K), '%% Date and time of this run: %s\n', ...
                params.ddmmyyyyhhmmss);
        fprintf(params.ccspectraAllSky_fid(K), '%s\n', ...
                '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(params.ccspectraAllSky_fid(K), '%s\t%s\t%s\t%s\t%s\n', ...
                '%start sec', 'theor sigma', 'freq (Hz)', ...
                'cc spec (real)', 'cc spec (imag)');
      end

    end

    if params.writeSensIntsToFiles
      % open files for sensitivity Integrand, etc. and write header info
      filename = [params.outputFilePrefix '_sensints.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.sensints_fid(K) = fopen(filename, 'w');
      
      fprintf(params.sensints_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.sensints_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.sensints_fid(K), '%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', ...
              'sens integrand (1/sec)');
      if params.doAllSkyComparison
        % open files for sensitivity Integrand, etc. and write header info
        filename = [params.outputFilePrefix '_AllSky_sensints.job' ...
                    num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
        params.sensintsAllSky_fid(K) = fopen(filename, 'w');
        
        fprintf(params.sensintsAllSky_fid(K), '%% Date and time of this run: %s\n', ...
                params.ddmmyyyyhhmmss);
        fprintf(params.sensintsAllSky_fid(K), '%s\n', ...
                '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(params.sensintsAllSky_fid(K), '%s\t%s\t%s\t%s\n', ...
                '%start sec', 'theor sigma', 'freq (Hz)', ...
                'sens integrand (1/sec)');
      end
    end
    
    if params.writeCoherenceToFiles
      % open files for the coherence, etc. and write header info
      filename = [params.outputFilePrefix '_coherence.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.coherence_fid(K) = fopen(filename, 'w');

      fprintf(params.coherence_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.coherence_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.coherence_fid(K), '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
              '%start sec', 'freq (Hz)', 'Coherence', ' real(csd) ',' imag(csd) ', '  psd1  ', '  psd2');
    end;
    
    if params.writeOptimalFiltersToFiles
      % open files for optimal filter and write header info
      filename = [params.outputFilePrefix '_optimal.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.optimal_fid(K) = fopen(filename, 'w');

      fprintf(params.optimal_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.optimal_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.optimal_fid(K), '%s\t%s\t%s\n', ...
              '%start sec', 'freq (Hz)', 'Optimal Filter Real / Imag');

      if params.doAllSkyComparison
        % open files for optimal filter and write header info
        filename = [params.outputFilePrefix '_AllSky_optimal.job' ...
                    num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
        params.optimalAllSky_fid(K) = fopen(filename, 'w');

        fprintf(params.optimalAllSky_fid(K), '%% Date and time of this run: %s\n', ...
                params.ddmmyyyyhhmmss);
        fprintf(params.optimalAllSky_fid(K), '%s\n', ...
                '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(params.optimalAllSky_fid(K), '%s\t%s\t%s\n', ...
                '%start sec', 'freq (Hz)', 'Optimal Filter Real / Imag');
      end

    end

    if params.writeOverlapReductionFunctionToFiles
      % open files for optimal filter and write header info
      filename = [params.outputFilePrefix '_orf.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.orf_fid(K) = fopen(filename, 'w');

      fprintf(params.orf_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.orf_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.orf_fid(K), '%s\t%s\t%s\n', ...
              '%start sec', 'freq (Hz)', 'OverlapReductionFunction Real / Imag');
    end

    if params.writeCalPSD1sToFiles
      % open files for PSD number 1 and write header info
      filename = [params.outputFilePrefix '_psd1.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.calpsd1_fid(K) = fopen(filename, 'w');

      fprintf(params.calpsd1_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.calpsd1_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.calpsd1_fid(K), '%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', 'calibrated PSD 1');
    end

    if params.writeCalPSD2sToFiles
      % open files for PSD number 2 and write header info
      filename = [params.outputFilePrefix '_psd2.job' ...
                  num2str(params.jobNumber) '.trial' num2str(K) '.dat'];
      params.calpsd2_fid(K) = fopen(filename, 'w');

      fprintf(params.calpsd2_fid(K), '%% Date and time of this run: %s\n', ...
              params.ddmmyyyyhhmmss);
      fprintf(params.calpsd2_fid(K), '%s\n', ...
              '% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      fprintf(params.calpsd2_fid(K), '%s\t%s\t%s\t%s\n', ...
              '%start sec', 'theor sigma', 'freq (Hz)', 'calibrated PSD 2');
    end

  end; % for K=1:params.numTrials

return;