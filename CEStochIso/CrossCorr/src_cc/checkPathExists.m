% A helper function to check for the existence of a path and
% standardise error messages.
%   desc  - a short description, generally just the name of the
%           variable to print out
%   p - the path to check (a string)
%
function checkPathExists(desc, p)

  if (~exist(p, 'dir'))
    error(sprintf('%s \''%s\'' not found.', desc, p));
  end;

return;
