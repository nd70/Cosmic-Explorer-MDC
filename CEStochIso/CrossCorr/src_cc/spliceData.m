function out = spliceData(which_buffer, x, y, N, which_IFO)

%  spliceData -- splices data contained in two arrays 
%
%  Splices data contained in two arrays (and possibly a buffer from a
%  previous call to the function) shifted relative to one another by 
%  half their length and weighted by a sinusoidal taper.
%
%  Inputs:
%    which_buffer: selects which buffer from an array of buffers
%    x: input array
%    y: input array (to be shifted to the right relative to x)
%    N: length of the input (and output) arrays
%    which_IFO: chooses between buffer1 and buffer2 (default=1)
%      (need two different sized buffers if the time-series for two
%       detectors have different lengths)
%
%  Outputs:
%    out: output array containing the spliced data
%
%  Persistent:
%    buffer1,2: input/output arrays containing data leftover from 
%      the previous call to spliceData (has half the length of x,y)
%
%  Routine written by Joseph D. Romano.
%
%  Contact Joseph.Romano@astro.cf.ac.uk
%
% $Id:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% persistent buffers for two IFOs
persistent buffer1
persistent buffer2
maxsize = 10;

% check that selected buffer is in range
if (which_buffer <= 0) | (which_buffer > maxsize) 
  error('selected buffer is out of range'); 
end

% check that length N >= 2
if N<2
  error('length N < 2')
end

% check that x, y have length N
if (length(x)~=N)
  error('length of x array does not equal N'); 
end
if (length(y)~=N)
  error('length of y array does not equal N'); 
end

% determine middle index (half of length) 
if mod(N,2)==0
  N_mid =  N/2;
else
  N_mid = (N-1)/2;
end

% determine which_IFO (default=1)
try
  which_IFO;
catch
  which_IFO = 1;
end

% check that selected IFO is in range
if (which_IFO~=1) & (which_IFO~=2)
  error('selected IFO is out of range');
end

% choose buffer appropriate for selected IFO, storing in temporary 
% buffer variable
if which_IFO==1
  buffer = buffer1;
else
  buffer = buffer2;
end

% initialize buffer to zeros if first pass through
if isempty(buffer)
  % first pass through; initialize buffer to 0's
  buffer = zeros(maxsize,N_mid);
end

% combine the two input arrays (offset by half their length) with one 
% another and with data leftover in the buffer 
if max(abs(buffer(which_buffer,1:N_mid)))==0
  % no data leftover in the buffer
  out(1:N_mid) = x(1:N_mid);  
else
  out(1:N_mid) = buffer(which_buffer,:) + ... 
                 x(1:N_mid).*sin(pi*[0:N_mid-1]/N);
end

out(N_mid+1:N) = x(N_mid+1:N).*sin(pi*[N_mid:N-1]/N) + ...
                 y(1:N-N_mid).*sin(pi*[0:N-N_mid-1]/N);

% save data for the next call
buffer(which_buffer,:) = y(N-N_mid+1:N).*sin(pi*[N-N_mid:N-1]/N);

% copy data from temporary buffer to persistent buffer for selected IFO
if which_IFO==1
  buffer1 = buffer;
else
  buffer2 = buffer;
end

% debugging plots
%figure(1); plot([1:N],x,'-*');
%figure(2); plot([1:N],y,'-*');
%figure(3); plot([1:N],out,'-*');
%figure(4); plot([1:N_mid],buffer(which_buffer,:),'-*');

% clear temporary buffer
clear buffer;

return
