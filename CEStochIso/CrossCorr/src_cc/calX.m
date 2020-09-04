function [X]=calX(C,P1,P2,H,glm,tsidereal,...
                  T,w1w2bar,w1w2squaredbar,mask,FreqIntFlag)

% calcluates the Fisher information matrix as
%
%                          H              w1w2bar^2
% X_lm' = T * gamma^*_lm ------- * C * ----------------
%                         P1 P2         w1w2squaredbar
%
% Assumes that C,P1,P2 and H are all of the same length and units (strain^2/Hz).
% In particular C has to be already corrected by 1/w1w2bar.
% This was NOT the case for old stochastic.m and C=s1^*s2,
%
% Input:
%
%   glm - a structure containing the following fields:
%
%   glm.data - values of gammaLM(f,tsidereal=0) for L=0,...,Lmax
%              and M=-L,..L for f = flow+deltaF*[0:numFreqs-1]
%              packed as a 2-d matrix as shown below:
% 
%              | flow .... fhigh
%   ---------------------------------
%   M=-L,L=Lmax|
%   .......    |
%   M=-1,L=Lmax|
%   .......    |
%   M=-1,L=1   |
%   M=0,L=0    |
%   M=0,L=1    |
%   .......    |
%   M=0,L=Lmax |
%   M=1,L=1    |
%   .......    |
%   M=1,L=Lmax |
%   .......    |
%   M=L=Lmax   |
%
%   glm.lvec     - array of l-values corresponding to the above packing
%   glm.mvec     - array of m-values corresponding to the above packing
%   glm.Lmax     - maximum value of L
%   glm.flow     - lowest frequency value (Hz)
%   glm.deltaF   - frequency spacing (Hz)
%   glm.numFreqs - number of discrete frequencies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract relevant quantities from glm structure
deltaF = glm.deltaF;
lvec = glm.lvec;
mvec = glm.mvec;

% construct phasor associated with sidereal time
phase1=i*pi*tsidereal/12;
phasor=exp(phase1*mvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  FreqIntFlag;
catch
  FreqIntFlag = 1;
end; 

if FreqIntFlag
   Xpos=...
  conj(glm.data) * (mask.*H./(P1.*P2) .* C) ...
  .* conj(phasor) ...
  .* (w1w2bar^2 / w1w2squaredbar *deltaF);
    
  % calculate X using the symmetry relation in l as well
  Xneg=...
   conj(glm.data) * (mask.*H./(conj(P1).*conj(P2)) .* conj(C)) ...
  .* conj(phasor) ...
  .* (w1w2bar^2 / w1w2squaredbar *deltaF) .*(-1).^lvec;
  X = T * (Xpos + Xneg);
else
   Xpos=...
  conj(glm.data) ...
  .* ((conj(phasor) .* (w1w2bar^2 / w1w2squaredbar *deltaF)) ...
      * transpose((mask.*H./(P1.*P2) .* C))) ;
    
  % calculate X using the symmetry relation in l as well
  Xneg=...
   conj(glm.data)  ...
  .* ((conj(phasor) .* (w1w2bar^2 / w1w2squaredbar *deltaF) .*(-1).^lvec) ...
     * transpose(mask.*H./(conj(P1).*conj(P2)) .* conj(C)));
  % Positive freqs from flow to fhigh, then negative freqs from -flow to -fhigh
  X = T * ([Xpos,Xneg]);
  %error('Frequency bin-by-bin not implemented yet');
end

