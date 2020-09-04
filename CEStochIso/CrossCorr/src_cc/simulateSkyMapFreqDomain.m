function [P1, P2, CP] = simulateSkyMapFreqDomain(Hf, intLog, ...
                                                 flow, deltaF, numFreqs, ...    
                                                 siderealtime, ...
                                                 glm, g1lm, g2lm, ...
                                                 det1, det2, ...
                                                 isSpH, coord1, coord2, map)
%
%  simulateSkyMap --- simulates a spatially-distributed unpolarized SGWB
%
%  arguments: 
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             flow       - low freq
%             deltaF     - delta F
%             numFreqs   - number of discrete frequencies
%             siderealtime - sidereal time in hours (0..24h)
%             glm, g1lm, g2lm - spherical harmonic components of the overlap reduction functions (for detectors 12, 11, 22)
%             det1,det2  - detector structures containing position r and tensor d
%             isSpH      - Type of map: true: map contains complex spherical harmonics; false: map is pixel map
%             coord1     - either vector of l or right ascension in hours
%                          takes vector for multiple point sources, must have same length as coord2
%             coord2     - either vector of m or declination in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as coord1
%             map        - map data; either complex spherical haromnics, or value at pixel
%
%  output:    P1,P2,CP   - power spectra and cross-power spectrum in the two detectors due to an SGWB
%
%  Routine written by Stefan Ballmer, Joe Romano.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert Hf into a frequency series appropriate for SpH2CrossPower and Map2CrossPower
H=checkArgumentFreqSeries(Hf,flow,deltaF,numFreqs,intLog);

% construct expected signal power spectra and cross-power
if isSpH == true
  % map is complex plms
  lvec = coord1;
  mvec = coord2;
  plm = map;

  % calculate expected signal power in an individual detector and 
  % the expected cross-power for the given plm's
  P1 = SpH2CrossPower(g1lm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  P1 = real(P1);
  P2 = SpH2CrossPower(g2lm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  P2 = real(P2);
  CP = SpH2CrossPower(glm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
else
  % pixel map
  ra = coord1;
  decl = coord2;
  RADECL = [ra decl];
  dOmega = 1;

  % calculate expected signal power in an individual detector and 
  % the expected cross-power for the given map
  P1 = Map2CrossPower(det1, det1, map, dOmega, RADECL, siderealTime, ...
                     numFreqs, flow, deltaF, H.data);
  P1 = real(P1);
  P2 = Map2CrossPower(det2, det2, map, dOmega, RADECL, siderealTime, ...
                     numFreqs, flow, deltaF, H.data);
  P2 = real(P2);
  CP = Map2CrossPower(det1, det2, map, dOmega, RADECL, siderealTime, ...
                     numFreqs, flow, deltaF, H.data);

end

% print out warning message if P1, P2, CP correspond to an unphysical signal
detC = P1.*P2-CP.*conj(CP);
if min(detC) <= 0
  warning('Trying to simulate an unphysical signal')
end

return
