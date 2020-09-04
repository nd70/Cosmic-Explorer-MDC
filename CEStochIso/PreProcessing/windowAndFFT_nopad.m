function y = windowAndFFT_nopad(x, window, fftLength)
%
%  windowAndFFT --- apply window and FFT
% 
%  windowAndFFT(x, window, fftLength) windows, 
%  does NOT zero-pads (to fftLength)
%  and DFTs a real time series of length N, sampling period deltaT, to
%  produce a complex frequency series of length fftLength/2.
% 
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%  modified by E. Thrane
%
%  $Id: windowAndFFT.m,v 1.6 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract metadata
N = length(x.data);
deltaT = x.deltaT;

% make sure that x and window have the same length
if ( length(x.data) ~= length(window) )
  error('size mismatch');
end

% make sure that fftLength >= length of the data 
if (fftLength < N)
  error('fftLength < length of data');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% zero-pad x
%eht: z = zeros(fftLength-N,1);
%eht: xbar = [ x.data.*window ; z ];
xbar = [x.data.*window];

% take the DFT of xbar 
y_temp = fft(xbar);

%% Behavior is different for real time series and for complex,
%% heterodyned time series

if ( isfield(x,'fbase') & (~ isnan(x.fbase)) )
  % put frequencies in ascending order
  offset = transpose(0:length(y_temp)-1) - length(y_temp);
  y_temp = fftshift(y_temp);
  offset = fftshift(offset);
  if ( isfield(x,'phase') & (~ isnan(x.phase)) )
    if ( x.phase ~= 0 )
      y_temp = y_temp .* exp(i*x.phase);
    end
  end
  y.symmetry = 0;
else
  % extract non-trivial information 
%eht:  y_temp = y_temp( 1:floor((fftLength/2 + 1)) );
  y_temp = y_temp( 1:floor((fftLength/2)) );
  y.symmetry = 1;
end

% normalize to include the dt which makes it an approximation to the CFT
y_temp = deltaT*y_temp;

% fill structure
y.data = y_temp;
y.deltaF = 1/(deltaT*fftLength);
if ( isfield(x,'fbase') & (~ isnan(x.fbase)) )
  y.flow = x.fbase + offset(1)*y.deltaF;
else
  y.flow = 0;
end

return

