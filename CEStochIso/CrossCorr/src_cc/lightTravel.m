function [deltaT, phi] = lightTravel(detectorPair)
%
% calculates the light travel time and angle between the separation 
% vector and the greenwich meridian for a pair of detectors.
%
% Input:
%
%   detectorPair - a string containing one of the following pairs of
%                  letters corresponding to different detector pairs 
%
%        HH: Hanford-Hanford
%        LL: Livingston-Livingston
%        VV: Virgo-Virgo
%        HL: Hanford-Livingston
%        HV: Hanford-Virgo
%        HG: Hanford-GEO
%        HT: Hanford-TAMA
%        LV: Livingston-Virgo
%        LG: Livingston-GEO
%        LT: Livingston-TAMA
%        VG: Virgo-GEO
%        VT: Virgo-TAMA
%        GT: GEO-TAMA
%
% Output:
%
%   deltaT - light travel time (in sec) between the two detectors
%   phi    - angle (in degrees) between the separation vector and 
%            greenwich meridian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure detgeom package is in path
%path(path,'../../../utilities/detgeom/matlab/');

% get site location and orientation information for each detector
switch detectorPair
  case 'HH'
    d1 = getdetector('LHO');
    d2 = getdetector('LHO');
  case 'LL'
    d1 = getdetector('LLO');
    d2 = getdetector('LLO');
  case 'VV'
    d1 = getdetector('VIRGO');
    d2 = getdetector('VIRGO');
  case 'HL'
    d1 = getdetector('LHO');
    d2 = getdetector('LLO');
  case 'HV'
    d1 = getdetector('LHO');
    d2 = getdetector('VIRGO');
  case 'HG'
    d1 = getdetector('LHO');
    d2 = getdetector('GEO600');
  case 'HT'
    d1 = getdetector('LHO');
    d2 = getdetector('TAMA');
  case 'LV'
    d1 = getdetector('LLO');
    d2 = getdetector('VIRGO');
  case 'LG'
    d1 = getdetector('LLO');
    d2 = getdetector('GEO600');
  case 'LT'
    d1 = getdetector('LLO');
    d2 = getdetector('TAMA');
  case 'VG'
    d1 = getdetector('VIRGO');
    d2 = getdetector('GEO600');
  case 'VT'
    d1 = getdetector('VIRGO');
    d2 = getdetector('TAMA');
  case 'GT'
    d1 = getdetector('GEO600');
    d2 = getdetector('TAMA');
  otherwise
    error('not a valid detector pair');
end

% calculate light travel time and angle between baseline and greenwich 
% meridian for the detector pair
speedOfLight = 299792458; % m/s
deltaX = d1.r - d2.r;
deltaT = norm(d1.r - d2.r)/speedOfLight;
if norm(deltaX)==0
  % co-located detectors
  phi = 0;
else
  phi = (180/pi)*atan(deltaX(2)/deltaX(1));
end

return

