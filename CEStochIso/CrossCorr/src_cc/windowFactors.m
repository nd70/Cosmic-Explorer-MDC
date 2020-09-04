function [w1w2bar, w1w2squaredbar, w1w2ovlsquaredbar] = ...
  windowFactors(window1, window2)
%
%  windowFactors -- calculate averages of products of windows
%
%  [w1w2bar, w1w2squaredbar, w1w2ovlsquaredbar] =
%  windowFactors(window1, window2) calculcates window factors needed 
%  to modify normalization factors, theoretical variances, error bars, 
%  etc. for non-trivial windowing.  
%
%  NOTE: The first two factors are needed even for non-overlapping
%  windows.
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: windowFactors.m,v 1.6 2007-03-21 16:07:19 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1 = length(window1);
N2 = length(window2);
Nred = gcd(N1,N2);

% check for incommensurate window lengths
if (Nred == 1)
  error('size mismatch\n');
end

% If window lengths are different, select reduced windows
% (points at corresponding times)
window1red = window1(1+(0:Nred-1)*N1/Nred);
window2red = window2(1+(0:Nred-1)*N2/Nred);

% extract 1st and 2nd half of windows
firsthalf1  = window1red(1:floor(Nred/2));
secondhalf1 = window1red(floor(Nred/2)+1:Nred);
firsthalf2  = window2red(1:floor(Nred/2));
secondhalf2 = window2red(floor(Nred/2)+1:Nred);
 
% calculate window factors
w1w2bar = mean(window1red.*window2red);
w1w2squaredbar = mean((window1red.^2) .* (window2red.^2));
w1w2ovlsquaredbar = ...
  mean((firsthalf1.*secondhalf1) .* (firsthalf2.*secondhalf2));

return

