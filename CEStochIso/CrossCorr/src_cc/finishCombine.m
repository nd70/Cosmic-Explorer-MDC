function finishCombine(params, K)

% finished the live combining of results of isotropic analysis
% and writes the values to files
% 
% input: params struct
% input: K - trial index
%
% Routine copied from stochastic.m and modified by Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: finishCombine.m,v 1.4 2007-08-09 17:46:44 vmandic Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  global COMB;

  %SHIFT: have to worry about different trials

  if params.doCombine

    if (params.doOverlap)
      COMB.ccVars_ovl(COMB.ovladdctr:end) = [];
      COMB.combined_ccstats =  COMB.combined_ccstats / sum(1./COMB.ccVars_ovl);
      COMB.combined_ccspectra =  COMB.combined_ccspectra / sum(1./COMB.ccVars_ovl);
      errorBar = sqrt(1/sum(1./COMB.ccVars_ovl));
    else % no overlap
      COMB.ccVars(COMB.addctr:end) = [];
      COMB.combined_ccstats =  COMB.add_ccstats / sum(1./COMB.ccVars);
      COMB.combined_ccspectra =  COMB.add_ccspectra / sum(1./COMB.ccVars);
      COMB.combined_sensint =  COMB.add_sensint;
      errorBar = sqrt(1/sum(1./COMB.ccVars)) / params.segmentDuration;
    end

    if (params.doBadGPSTimes)
      COMB.badGPSTimes(COMB.badGPSTimesctr:end) = [];
    end

    % take out the scaling by segment duration
    COMB.combined_ccstats = COMB.combined_ccstats * params.segmentDuration;
    errorBar = errorBar * params.segmentDuration;
    COMB.combined_ccspectra = COMB.combined_ccspectra * params.segmentDuration;
    COMB.combined_sensint = COMB.combined_sensint / params.segmentDuration^2;

    % Write results to files
    % TODO: Seems to be always done if doCombine is chosen, perhaps there should be a user option
    if (length(COMB.combined_ccspectra) == params.numFreqs) && ...
          (length(COMB.combined_sensint  ) == params.numFreqs)
      
      if (params.writeOutputToMatFile)
        % For matfiles, we copy the important parts of the COMB data structure, keeping exactly
        % the same format of the output structure as we do for results for each segment so that
        % post-processing can be done in the same way. The combined data is written to a
        % separate matfile with the same name apart from having '_combined' inserted before
        % the job number.
        % We treat the startTime as the segmentStartTime for compatability with post-processing
        params.combMatFile.segmentStartTime(1, 1) = params.startTime;
        % errorBar is read as ccSigma in post-processing
        params.combMatFile.ccSigma(1, K) = errorBar;
        params.combMatFile.ccStat(1, K) = COMB.combined_ccstats;
        params.combMatFile.ccSpec(1, K) = constructFreqSeries(COMB.combined_ccspectra, params.flow, params.deltaF);
        params.combMatFile.sensInt(1, K) = constructFreqSeries(COMB.combined_sensint, params.flow, params.deltaF);
        if params.doBadGPSTimes
          % We can't write the bad times array directly to the matfile so we need to read in the cellarray
          % of bad GPS times and then write it out again.
          % The times are listed as a column vector indexed by trial in a cellarray badGPSTimes{K}
          badGPSTimes = params.combMatFile.badGPSTimes;
          badGPSTimes{K} = COMB.badGPSTimes.';
          params.combMatFile.badGPSTimes = badGPSTimes;
        end;
      else % write to text file
        fprintf(params.combined_ccstats_fid(K),'%9.1f\t%e\t%e\n', params.startTime, ...
                COMB.combined_ccstats, errorBar);
        
        fprintf(params.combined_ccspectra_fid(K), '%9.1f\t%e\t%e\t%e\t%e\n', ...
                [ params.startTime*transpose(ones(params.numFreqs, 1)); ...
                  errorBar*transpose(ones(params.numFreqs, 1)); transpose(params.f); ...
                  transpose(real(COMB.combined_ccspectra)); ...
                  transpose(imag(COMB.combined_ccspectra)) ]);
        
        fprintf(params.combined_sensints_fid(K), '%9.1f\t%e\t%e\t%e\n', ...
                [ params.startTime*transpose(ones(params.numFreqs, 1)); ...
                  errorBar*transpose(ones(params.numFreqs, 1)); ...
                  transpose(params.f); transpose(COMB.combined_sensint) ]);

        % For text output we write the badGPSTimes to a separate file, but for
        % matfile output they are included in the .mat file
        if params.doBadGPSTimes
          fprintf(params.badGPSTimes_fid(K),'%9.1f\n', transpose(COMB.badGPSTimes));
        end;
      end; % if params.writeOutputToMatFile

    else % lengths of spectra and numFreqs don't match for some reason
      warning('No valid ccspectra or sensint data, combined_ccspectra and combined_sensints not saved.');
    end; % if writing combined data to files
  
  end % if params.doCombine
  
return;