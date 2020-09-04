function [Gamma]=calFisher(P1,P2,H,glm,tsidereal,...
                           T,w1w2bar,w1w2squaredbar,mask,FreqIntFlag)

% calcluates the Fisher information matrix as
%
%                                  H^2                  w1w2bar^2
% Gamma_lm_l'm' = T * gamma^*_lm ------- gamma_l'm' * --------------
%                                 P1 P2               w1w2squaredbar
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Fisher matrix is never calculated bin-by-bin
%if FreqIntFlag
if true
  Gamma=...
  conj(glm.data) * ( ( (mask.*H.^2./(P1.*P2))*ones(1,size(glm.data,1)) ) .* transpose(glm.data) ) ...
  .* (conj(phasor) * transpose(phasor)) ...
  .* (w1w2bar^2 / w1w2squaredbar *deltaF);
else
  error('Frequency bin-by-bin not implemented yet');
end

% to sum over negative frequencies as well, exploiting symmetries of
% gammaLM
negfreqfactor=((-1).^lvec)*((-1).^transpose(lvec));
Gamma = T * Gamma .* (ones(length(lvec),length(lvec)) + negfreqfactor);
