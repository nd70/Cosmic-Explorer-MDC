function [y, index1, index2, frac1, frac2]=coarseGrain(x, flowy, deltaFy, Ny)
%
%  coarseGrain --- coarse grain a frequency-series 
%
%  coarseGrain(x, flowy, deltaFy, Ny) returns a frequency-series structure 
%  coarse-grained to the frequency values f = flowy + deltaFy*[0:Ny-1].  
%
%  coarseGrain also returns the indices of the lowest and highest frequency 
%  bins of x that overlap with y (0 <= index1 <= length(x); 
%  1 <= index2 <= length(x)+1) and the fractional contributions from these 
%  frequency bins.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: coarseGrain.m,v 1.4 2005-02-24 14:36:27 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract metadata
Nx = length(x.data);
flowx = x.flow;
deltaFx = x.deltaF;
% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that frequency spacings and number of elements are > 0
if ( (deltaFx <= 0) | (deltaFy <= 0) | (Ny <= 0) )
  error('bad input arguments');
end

% check that freq resolution of x is finer than desired freq resolution
if ( deltaFy < deltaFx )
  error('deltaF coarse-grain < deltaF fine-grain');
end

% check desired start frequency for coarse-grained series
if ( (flowy - 0.5*deltaFy) < (flowx - 0.5*deltaFx) )
  error('desired coarse-grained start frequency is too low');
end

% check desired stop frequency for coarse-grained series
fhighx = flowx + (Nx-1)*deltaFx;
fhighy = flowy + (Ny-1)*deltaFy;

if ( (fhighy + 0.5*deltaFy) > (fhighx + 0.5*deltaFx) )
  error('desired coarse-grained stop frequency is too high');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indices for coarse-grained series
i = transpose([0:1:Ny-1]);

% calculate the low and high indices of fine-grained series 
jlow  = 1 + floor1( (flowy + (i-0.5)*deltaFy - flowx - 0.5*deltaFx)/deltaFx ); 
jhigh = 1 + floor1( (flowy + (i+0.5)*deltaFy - flowx - 0.5*deltaFx)/deltaFx ); 

index1 = jlow(1);
index2 = jhigh(end);

% calculate fractional contributions
fraclow  = (flowx + (jlow+0.5)*deltaFx - flowy - (i-0.5)*deltaFy)/deltaFx;
frachigh = (flowy + (i+0.5)*deltaFy - flowx - (jhigh-0.5)*deltaFx)/deltaFx;
  
frac1 = fraclow(1);
frac2 = frachigh(end);

% calculate coarse-grained values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum of middle terms
% NOTE:  since this operation takes more time than anything else, i made 
% it into separate routine called sumTerms() which is then converted into
% a MEX file
jtemp = jlow+2;
midsum = sumTerms(x.data, jtemp, jhigh);

% calculate all but final value of y
ya = (deltaFx/deltaFy)*(x.data(jlow (1:Ny-1)+1).*fraclow (1:Ny-1) + ...
                        x.data(jhigh(1:Ny-1)+1).*frachigh(1:Ny-1) + ...
                        midsum(1:Ny-1) );

% calculate final value of y
if (jhigh(Ny) > Nx-1)
  % special case when jhigh exceeds maximum allowed value
  yb = (deltaFx/deltaFy)*(x.data(jlow(Ny)+1)*fraclow(Ny) + midsum(Ny));

else
  yb = (deltaFx/deltaFy)*(x.data(jlow (Ny)+1)*fraclow (Ny) + ...
                          x.data(jhigh(Ny)+1)*frachigh(Ny) + ...
                          midsum(Ny) );
end

% fill structure for coarsed-grained frequency series
y.data = [ya; yb];
y.flow = flowy;
y.deltaF = deltaFy;
try
  y.symmetry = x.symmetry;

 
 catch 
 

  y.symmetry = 0;
end

return

function indexC = floor1(index)
% function floor1(index)
% modified floor function to account for precision error 
  if abs(index-round(index)) < 1e-8
    indexC = round(index);
  else
    indexC = floor(index);
  end
return
