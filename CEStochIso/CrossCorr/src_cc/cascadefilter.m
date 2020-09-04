function y = cascadefilter(x, n, fknee, fr)
% function y = cascadefilter(x, n, fknee, fr)
% This function applies an order-n Butterworth high-pass filter by applying
% n/2 consecutive second-order Butterworth filters.  Suggestion from Ed Daw.
%   x = input data
%   n = order
%   fknee = knee frequency
%   fr (re)sample rate
%   y = ouput data
% Eric Thrane

% only run with numbers of n divisible by 4
if mod(n,4)~=0
  error('n must be divisible by 4');
end

% get n=4 filter coefficients
[b,a] = butter(4, fknee/(fr/2), 'high');

% number of n=4 filters to apply
m = n/4;

% apply filters
for ii=1:m
  x = filtfilt(b, a, x);
end

% finish up
y = x;

return
