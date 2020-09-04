function params=closeOutputFiles(params)
%
% Close the output files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (params.writeOutputToMatFile)
    % There is no 'close' function for matfiles but we can make it unwritable
    params.matFile.Properties.Writable = false;
  else
    params = closeOutputTextFiles(params);
  end;

return;
