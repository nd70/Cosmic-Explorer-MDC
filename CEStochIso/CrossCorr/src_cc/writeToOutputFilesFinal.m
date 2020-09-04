function writeToOutputFilesFinal(params, K)
%
% Adds any final output to the output files that doesn't depend on the segment index.
% 
% input: params struct
%        K - trial number
%
% $Id: writeToOutputFiles.m,v 1.2 2007-07-18 01:18:34 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (params.writeOutputToMatFile)
    
    % coh
    % There is only one coherence for each run so we don't use the segment index,
    % but there will be multiple trials
    if params.writeCoherenceToFiles
      params.matFile.coh(1, K) = params.coh;
    end

    % gamma
    if params.writeOverlapReductionFunctionToFiles
      % write overlap reduction function to output files
      params.matFile.gamma(1, K) = params.gamma;
    end

  else % write to text files

    if (params.writeCoherenceToFiles)
      % write coherence to files (only write once per job)
      fprintf(params.coherence_fid(K), '%d\t%e\t%e\t%e\t%e\t%e\t%e\n', ...
              [ params.coh.coherenceStartTime*ones(1, length(params.coh.f_coh)); ...
                params.coh.f_coh; params.coh.coh_avg'; real(params.coh.cavg'); ...
                imag(params.coh.cavg'); params.coh.p1avg'; params.coh.p2avg' ]);
    end;

  end; % if writeOutputToMatFile

return;