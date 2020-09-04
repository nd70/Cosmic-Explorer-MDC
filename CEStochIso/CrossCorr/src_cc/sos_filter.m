% function [x2] = sos_filter(x, n, fknee, fs)
%
% This function applies an order-n sos high-pass filter by applying second-order Butterworth filters.
%   x = input data
%   n = order
%   fknee = knee frequency
%   fs = sampling rate
%   x2 = output data

function x2=sos_filter(x, n, fknee, fs)

  [z, p, k] = butter(n, fknee/(fs/2), 'high');
  [sos, G] = zp2sos(z, p, k);
  x2 = filtfilt(sos, G, x);

return;
