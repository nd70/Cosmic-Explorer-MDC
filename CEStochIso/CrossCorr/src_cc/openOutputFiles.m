function params=openOutputFiles(params)

% Opens output files for writing
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: openOutputFiles.m,v 1.1 2007-06-26 00:18:48 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (params.writeOutputToMatFile)
    params = openOutputMatFile(params);
  else
    params = openOutputTextFiles(params);
  end;

end
