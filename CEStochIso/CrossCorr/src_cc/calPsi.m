function [Psi]=calPsi(P1,P2,H,glm,tsidereal,...
                      T,w1w2bar,w1w2squaredbar,mask,FreqIntFlag)

% calcluates the correction factor to the Fisher matrix for 
% simulatedSkyMap signal injections
%
%                                                             H^2 H_inj^2                 w1w2bar^2
% Psi_lm_l'm' = T int_df (-1)^l' gamma^*_lm (pdotgamma_inj)^2 ----------- gamma_l'm' * --------------
%                                                              P1^2 P2^2                w1w2squaredbar
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global SKY_MAP_MEMORY;

% extract relevant quantities from glm structure (refers to filter)
lvec = glm.lvec;
mvec = glm.mvec;
flow = glm.flow;
deltaF = glm.deltaF;
numFreqs = glm.numFreqs;

% construct phasor associated with sidereal time
phase1=i*pi*tsidereal/12;
phasor=exp(phase1*mvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  FreqIntFlag;
catch
  FreqIntFlag = 1;
end; 

% extract plm's, gammalm's, and power spectrum of injected sky map
plm_inj = SKY_MAP_MEMORY.map;
glm_inj = SKY_MAP_MEMORY.glm;
Hf = checkArgumentFreqSeries(SKY_MAP_MEMORY.Hf,flow,deltaF,numFreqs,SKY_MAP_MEMORY.intLog);
H_inj = Hf.data;

% calculate dot product of injected plm's and gamma_lm's wrt lm indices
pdotgamma_inj = transpose( transpose(plm_inj) * glm_inj.data );

% interpolate dot product to agree with filter discrete frequencies  
f_inj = glm_inj.flow + glm_inj.deltaF*[0:glm_inj.numFreqs-1]';
f = flow + deltaF*[0:numFreqs-1]';
pdotgamma_inj = exp(interp1(log(f_inj), log(pdotgamma_inj), log(f)));

% calculate gamma_lm's for negative frequency
glm_neg = zeros(size(glm.data,1),size(glm.data,2));
for kk=1:size(glm.data,2)
  glm_neg(:,kk) = ((-1).^lvec) .* glm.data(:,kk);
end

% calculate Psi
if FreqIntFlag
  Psi=...
  conj(glm.data) * ( ( (mask.* (pdotgamma_inj.^2) .* (H_inj.^2) .* H.^2./(P1.*P2).^2)*ones(1,size(glm.data,1)) ) .* transpose(glm_neg) ) ...
  .* (conj(phasor) * transpose(phasor)) ...
  .* (w1w2bar^2 / w1w2squaredbar *deltaF);
else
  error('Frequency bin-by-bin not implemented yet');
end

% to sum over negative frequencies as well, exploiting symmetries of gamma_lm
negfreqfactor=((-1).^lvec)*((-1).^transpose(lvec));
Psi = T * Psi .* (ones(length(lvec),length(lvec)) + negfreqfactor);
