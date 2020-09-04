% A helper function to check for the existence of a file and
% standardise error messages.
%   desc  - a short description, generally just the name of the
%           variable to print out
%   p - the full path to the file (a string)
%
function checkFileExists(desc, p)

  % First make sure we don't have a directory, which is also an error
  % (matlab treats directories as files too so just checking for
  % existence of a file is not enough)
  if (exist(p, 'dir'))
    error(sprintf('%s \''%s\'' is a directory.', desc, p));
  end;

  if (~exist(p, 'file'))
    error(sprintf('%s \''%s\'' not found.', desc, p));
  end;

return;