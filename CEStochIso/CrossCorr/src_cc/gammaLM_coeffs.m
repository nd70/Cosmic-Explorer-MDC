function gammaLM_coeffs(detectorPair, Lmax)
%
% Produces a .mat file containing the spherical bessel function 
% coefficients for the gammaLMs for a pair of detectors as calculated 
% by the Mathematica notebook gammaLM_coeffs.nb.  
%
% Input:
%   detectorPair - a string containing one of the following pairs of
%     letters corresponding to different detector pairs 
%
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
%        HH: Hanford-Hanford
%        LL: Livingston-Livingston
%        VV: Virgo-Virgo
%
%   Files of the form 
%
%      XXjcoeffsLMreal.dat, XXjcoeffsLMimag.dat   
%
%   where XX is one of the detector pairs above, L=0, 1, ..., Lmax, 
%   M=0, 1, ..., L.  These files contain the spherical bessel function 
%   coefficients for a particular gammaLM.  These files were originally 
%   created by the Mathematica notebook gammaLM_coeffs.nb.
%
% Output:
%   .mat file of the form
%
%      XXgammaLM_coeffs.mat
%
%   where XX is one of the detector pairs above.  This file contains 
%   the the following 2-d matrix of values
% 
%   L=0    M=0     coeff1   coeff2   coeff3  ...
%   L=1    M=0     coeff1   coeff2   coeff3  ...
%   L=1    M=1     coeff1   coeff2   coeff3  ...
%   ...
%
%   L=Lmax M=Lmax  coeff1   coeff2   coeff3  ...
% 
% Note: The spherical bessel function coefficients and the definition 
% of the gammaLM's as calculated by the Mathematica notebook follow the 
% conventions and definitions in Allen-Ottewill, PRD 56, 545 1997.  
% Renormalisation of the gammaLM's in both amplitude and phase may be 
% needed for subsequent analyses.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some common quantities
numCoeffsMax =  2 + floor((Lmax + 2)/2);
nRows = (Lmax+1)*(Lmax+2)/2; 
nCols = numCoeffsMax + 2;

% initialise data matrix
coeffs = zeros(nRows, nCols);

for L=0:Lmax

  for M=0:L

    % read in coeffs from files produced by mathematica notebook
    filename_real = [detectorPair 'jcoeffs' num2str(L) num2str(M) 'real.dat'];
    coeffs_real = ...
     textread(filename_real, '%f\n', -1, 'commentstyle', 'matlab');

    filename_imag = [detectorPair 'jcoeffs' num2str(L) num2str(M) 'imag.dat'];
    coeffs_imag = ...
     textread(filename_imag, '%f\n', -1, 'commentstyle', 'matlab');

    % fill-in matrix of coefficients
    numCoeffs = 2 + floor((L+2)/2);
    row = L*(L+1)/2 + 1 + M; 

    % coeffs for gammaLM
    coeffs(row,1) = L;
    coeffs(row,2) = M;
    coeffs(row,3:3+(numCoeffs-1)) = ...
        transpose(coeffs_real+sqrt(-1)*coeffs_imag);

  end

end

% save data array to a .mat file
fName = [detectorPair 'gammaLM_coeffs.mat'];
vName = 'coeffs';
save(fName,vName);

return
