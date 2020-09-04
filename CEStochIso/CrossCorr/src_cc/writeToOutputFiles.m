function writeToOutputFiles(params, I, K, segmentStartTime, result)
%
% Adds the results for the current segment and trial to the output files
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: writeToOutputFiles.m,v 1.2 2007-07-18 01:18:34 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Display results on screen if desired
  if params.writeResultsToScreen

    if params.doOverlap
      segmentNumber = (I-1) + (2*params.midSegment-1);
    else
      segmentNumber = (I-1) + params.midSegment; 
    end

    % display cc stats, theoretical sigmas to screen
    if (params.doDirectional)
      for k = 1:size(result.ccStat.data)
        fprintf('segment=%d, trial =%d, start sec = %d, CC stat = %e, theor sigma = %e\n', segmentNumber, K, ...
                segmentStartTime, result.ccStat.data(k), sqrt(result.ccVar));
      end;
    else
      if isreal(result.ccStat)
        fprintf('segment=%d, trial =%d, start sec = %d, CC stat = %e, theor sigma = %e\n', segmentNumber, K, ...
                segmentStartTime, result.ccStat, sqrt(result.ccVar));
      elseif (imag(result.ccStat)>=0)
        fprintf('segment=%d, trial =%d, start sec = %d, CC stat = %e + %e i, theor sigma = %e\n', segmentNumber, K, ...
                segmentStartTime, real(result.ccStat), imag(result.ccStat), sqrt(result.ccVar));
      else
        fprintf('segment=%d, trial =%d, start sec = %d, CC stat = %e - %e i, theor sigma = %e\n', segmentNumber, K, ...
                segmentStartTime, real(result.ccStat), -imag(result.ccStat), sqrt(result.ccVar));
      end;
    end;

  end; % params.writeResultsToScreen

  if params.writeOutputToMatFile
    writeOutputToMatFile(params, I, K, segmentStartTime, result);
  else
    writeOutputToTextFiles(params, I, K, segmentStartTime, result);
  end;
  
return;