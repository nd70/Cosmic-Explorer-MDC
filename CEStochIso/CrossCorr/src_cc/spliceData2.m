function out = spliceData(which_buffer, x, y, N)

%  spliceData -- splices data contained in two arrays 
%
%  Splices data contained in two arrays (and possibly a buffer from a
%  previous call to the function) shifted relative to one another by 
%  half their length and weighted by a sinusoidal taper.
%
%  Inputs:
%    which_buffer
%    x: input array
%    y: input array (to be shifted to the right relative to x)
%    N: length of the input (and output) arrays
%
%  Outputs:
%    out: output array containing the spliced data
%
%  Persistent:
%    buffer2: input/output array containing data leftover from the
%      previous call to spliceData (has half the length of x,y)
%
%  Routine written by Joseph D. Romano.
%
%  Contact Joseph.Romano@astro.cf.ac.uk
%
% $Id:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% persistent buffer
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

% initialize buffer to zeros
if isempty(buffer2)
  % first pass through; initialize buffer to 0's
  buffer2 = zeros(maxsize,N_mid);
end

% combine the two input arrays (offset by half their length) with one 
% another and with data leftover in the buffer 
if max(abs(buffer2(which_buffer,1:N_mid)))==0
  % no data leftover in the buffer
  out(1:N_mid) = x(1:N_mid);  
else
  out(1:N_mid) = buffer2(which_buffer,:) + ... 
                 x(1:N_mid).*sin(pi*[0:N_mid-1]/N);
end

out(N_mid+1:N) = x(N_mid+1:N).*sin(pi*[N_mid:N-1]/N) + ...
                 y(1:N-N_mid).*sin(pi*[0:N-N_mid-1]/N);

% save data for the next call
buffer2(which_buffer,:) = y(N-N_mid+1:N).*sin(pi*[N-N_mid:N-1]/N);

% debugging plots
%figure(1); plot([1:N],x,'-*');
%figure(2); plot([1:N],y,'-*');
%figure(3); plot([1:N],out,'-*');
%figure(4); plot([1:N_mid],buffer2(which_buffer,:),'-*');

return
