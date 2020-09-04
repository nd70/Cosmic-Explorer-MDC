function plotGammaLM(detectorPair, L, M, flow, fhigh)
%
% Plots gammaLM(f,t=0) for a given detector pair and for a particular 
% value of L,M using conventions and definitions of the gammaLM's from 
% Allen-Ottewill, PRD 56, 545 1997.  In particular, t=0 corresponds
% to 38.2 degrees east of Greenwich Meridian for Hanford-Livingston.
%
% Input:
%   detectorPair - a string containing one of the following pairs of
%     letters corresponding to different detector pairs 
%
%        HH: Hanford-Hanford
%        LL: Livingston-Livingston
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
%   flow, fhigh - min and max frequencies for gammaLM(f)
%   XXgammaLM_coeffs.mat - where XX is one of the detector pairs above
%     is a file containing the spherical bessel function
%     coefficients for the gammaLM following Allen-Ottewill, PRD 56,
%     545 1997 conventions.  The coefficients are arranged in a
%     matrix as follows:
%
%   coeffs: 
%     L=0    M=0     coeff1   coeff2   coeff3  ...
%     L=1    M=0     coeff1   coeff2   coeff3  ...
%     L=1    M=1     coeff1   coeff2   coeff3  ...
%     ...
%
%     L=Lmax M=Lmax  coeff1   coeff2   coeff3  ...
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','MATLAB:divideByZero');

% check that values of L,M are valid
if abs(M)>L
  error('Absolute value of M > L');
end

% plotting defaults
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontWeight','bold');
%set(0,'DefaultTextFontSize',12);
%set(0,'DefaultTextFontWeight','normal');
%set(0,'DefaultUicontrolFontSize',12);
%set(0,'DefaultUicontrolFontWeight','normal');
%set(0,'DefaultAxesTickLength',[0.02 0.01]);
set(0,'DefaultAxesLineWidth', 2.0);
set(0,'DefaultLineLineWidth', 2.0);
set(0,'DefaultAxesGridLineStyle','--');
set(0,'DefaultAxesXColor',[0 0 0]);
set(0,'DefaultAxesYColor',[0 0 0]);

% frequencies
numFreqs = 1000;
DeltaF = (fhigh-flow)/numFreqs;
f = transpose(flow + DeltaF*[0:numFreqs-1]);

% calculate light travel time between the detectors
[DeltaT, phi] = lightTravel(detectorPair);

% change of variables
x = 2*pi*f*DeltaT;
% special case for coincident-coaligned detectors
if DeltaT==0
  % need non-zero to get limiting value for jn(x)/x^n
  x=(1e-9)*ones(size(f)); 
end

% load gammaLM coefficients from appropriate .mat file
fName = [detectorPair 'gammaLM_coeffs.mat'];
vName = 'coeffs';
load(fName,vName);
if L > (2*size(coeffs,1)+0.25)^0.5 -1.5
  error('Requested value of L is too large');
end

% extract coefficients for M>0
numCoeffs =  2 + floor((L + 2)/2);
row = L*(L+1)/2 + 1 + abs(M);
jcoeffs = transpose(coeffs(row,3:3+(numCoeffs-1)));

% construct gammaLM for M>0
gammaLM = zeros(numFreqs,1);
for n=1:numCoeffs
  gammaLM = gammaLM + jcoeffs(n)*sphericalbessel(n-1+mod(L,2),x)./x.^(n-1);
end

% M<0
if M<0
  gammaLM= ((-1)^(L+abs(M))) * conj(gammaLM);
end
 
% make plot
figure(L*(L+1)/2+abs(M)+1);
plot(f, real(gammaLM), f, imag(gammaLM));
axis([flow, fhigh, -.3, .3]);
grid on;
xlabel('freq (Hz)');
ylabel('gammaLM(f)');
titleString = ['gammaLM(f,t=0) for ' detectorPair];
title(titleString);
legend(['re, l=',num2str(L),', m=',num2str(M)],...
       ['im, l=',num2str(L),', m=',num2str(M)]);
  
%feps=[detectorPair 'plotGamma' num2str(L) num2str(M) '.eps'];
%print('-depsc2',feps);

return

