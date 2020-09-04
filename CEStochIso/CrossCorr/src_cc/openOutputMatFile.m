function params=openOutputMatFile(params)

% Opens output file for writing in Matlab format. When writing output
% in matfile format the output is all stored in a single file. Since
% the contents of the file are updated 'in place' any existing files
% must be deleted at the start of the run to maintain consistent data
% within the file.
% 
% input and output: params struct
%
% $Id: openOutputFiles.m,v 1.1 2007-06-26 00:18:48 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Mat-file for all results ie. not combined
  matFileName = [ params.outputFilePrefix '.job' num2str(params.jobNumber) '.mat' ];

  % We need to delete the old output file if it exists
  if (exist(matFileName, 'file'))
    delete(matFileName);
  end;

  % Create an empty file and store the handle in params. The name is not
  % stored directly but can be recovered from params.matFile.Properties.Source
  params.matFile = matfile(matFileName, 'Writable', true);

  % Write current params for later reference. We first need to remove the large data
  % arrays that would otherwise bloat the save file
  fields_to_remove = { 'detector1', 'detector2', 'filt1', 'filt2', 'psd1', 'psd2', ...
                      'fft1', 'fft2', 'mask', 'f', 'gamma', 'cal1', 'cal2', 'matFile' };
  params.matFile.params = rmfield(params, fields_to_remove);

  % Insert an empty badGPSTimes cell array. Any bad times are stored here as a
  % cell array of vectors, where eg. badGPSTimes{K} is the array of bad times
  % in trial K
  if (params.doBadGPSTimes)
    params.matFile.badGPSTimes = cell(1, params.numTrials);
  end;

  % Mat-file for all combined results
  if (params.doCombine)
    combMatFileName = [ params.outputFilePrefix '_combined.job' num2str(params.jobNumber) '.mat' ];

    % We need to delete the old combined output file if it exists
    if (exist(combMatFileName, 'file'))
      delete(combMatFileName);
    end;

    params.combMatFile = matfile(combMatFileName, 'Writable', true);

    % Store params in combined data files as well
    params.combMatFile.params = rmfield(params, fields_to_remove);

    % Insert an empty badGPSTimes array
    if (params.doBadGPSTimes)
      params.combMatFile.badGPSTimes = cell(1, params.numTrials);
    end;

  end;

return;
