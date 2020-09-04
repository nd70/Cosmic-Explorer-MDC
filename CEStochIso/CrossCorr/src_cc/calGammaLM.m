function glm = calGammaLM(gammaLM_coeffsPath, detectorPair, Lmax, numFreqs, flow, deltaF)

% Calculates the gammaLM(f,t=0) for Greenwich mean sidereal time=0,
% for a detector pair using the conventions and definitions of the 
% gammaLM's from LIGO-T070045-00-U. This differs from Allen-Ottewill, 
% PRD 56, 545 1997 by an overall normalisation and phase factor.
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
%   Lmax     - max value of L
%   flow     - lowest frequency value (Hz)
%   deltaF   - frequency spacing (Hz)
%   numFreqs - number of discrete frequencies 
%
%   XXgammaLM_coeffs.mat - where XX is one of the detector pairs above
%     is a file containing the spherical bessel function coefficients 
%     for the gammaLM, following the 
%     conventions and definitions in Allen-Ottewill, PRD 56, 545 1997.
%
%     The coefficients are packed in a 2-d matrix having the form 
%     shown below:
%
%     L=0    M=0     coeff1   coeff2   coeff3  ...
%     L=1    M=0     coeff1   coeff2   coeff3  ...
%     L=1    M=1     coeff1   coeff2   coeff3  ...
%     ...
%
%     L=Lmax M=Lmax  coeff1   coeff2   coeff3  ...
% 
%
% Output:
%
%   A structure glm containing the following fields:
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
warning('off','MATLAB:divideByZero');

% load gammaLM coefficients from appropriate .mat file
if length(gammaLM_coeffsPath)==0
  fName = [detectorPair,'gammaLM_coeffs.mat'];
else
  fName = [gammaLM_coeffsPath,'/',detectorPair,'gammaLM_coeffs.mat'];
end
vName = 'coeffs';
load(fName,vName);
if Lmax > (2*size(coeffs,1)+0.25)^0.5 -1.5
  error('Requested value of Lmax is too large');
end

% calculate the light travel and the angle between the separation
% vector connecting the two detectors and the greenwich meridian
[DeltaT, phi] = lightTravel(detectorPair);

% frequency array
f = flow + deltaF*[0:numFreqs-1];
fhigh = f(end);

% change of variables
x = 2*pi*f*DeltaT;
% special case for coincident-coaligned detectors
if DeltaT==0
  % need non-zero to get limiting value for jn(x)/x^n
  x=(1e-9)*ones(size(f)); 
end

% pre-calculate functions that multiply the spherical bessel coeffs
% to improve efficiency
numCoeffsMax = 2 + floor((Lmax + 2)/2);
j_matrix = zeros(numCoeffsMax+1,numFreqs);
x_matrix = zeros(numCoeffsMax+1,numFreqs);
for n=1:numCoeffsMax+1
  j_matrix(n,:) = sphericalbessel(n-1,x);
  x_matrix(n,:) = x.^(n-1);
end
jbyx_e = j_matrix(1:numCoeffsMax,:)./x_matrix(1:numCoeffsMax,:);
jbyx_o = j_matrix(2:numCoeffsMax+1,:)./x_matrix(1:numCoeffsMax,:);

% initialize variables for positive values of m
numRows = (Lmax+1)*(Lmax+2)/2;
lvec = zeros(numRows,1);
mvec = zeros(numRows,1);
gammaLM = zeros(numRows,numFreqs);

for M=0:Lmax

  %fprintf('Working on M=%d out of %d\n',M,Lmax);

  % index into the gammaLM matrix for the M=L values (M=0,1, ..., Lmax)
  i0 = M*(Lmax+1)- M*(M-1)/2 + 1;

  % vector of m values, l values
  mvec(i0:(i0+Lmax-M)) = M;
  lvec(i0:(i0+Lmax-M)) = M:Lmax;
  for L=M:Lmax 

    % row index into the gammaLM matrix
    ii = i0 + L-M;

    % extract spherical bessel function coefficients for the gammaLMs
    numCoeffs =  2 + floor((L + 2)/2);
    row = L*(L+1)/2 + 1 + M;
    jcoeffs = coeffs(row,3:3+(numCoeffs-1));

    % construct gammaLM
    if mod(L,2)==0
      gammaLM(ii,:) = jcoeffs*jbyx_e(1:numCoeffs,:);
    else
      gammaLM(ii,:) = jcoeffs*jbyx_o(1:numCoeffs,:);
    end

    % rotate so that zero sidereal time corresponds to the greenwich
    % meridian and normalise to get agreement with LIGO-T070045-00-U
    phasefactor = exp(sqrt(-1)*M*phi*pi/180);
    normfactor = 4*pi/5;
    % NEED TO CHECK IF THESE ARE CORRECT NORMALISATION AND PHASE!!
    % Verified both norm and phase factor
    % Testfunction: testSignalFromNorthpole.m
    % Stefan Ballmer 11/19/2007
    gammaLM(ii,:) = gammaLM(ii,:)*normfactor*phasefactor;

  end

end

% extend gammaLM matrix to have values for negative m
gammaLMM = conj(gammaLM).*((-1).^(lvec+mvec) * ones(1,size(gammaLM,2)));
gammaLM = [flipud(gammaLMM(2+Lmax:end,:)); gammaLM]; 

% extend vectors of l and m values to negative values of m
lvec = [ flipud(lvec(2+Lmax:end)); lvec];
mvec = [-flipud(mvec(2+Lmax:end)); mvec];

% assign values to the output glm structure
glm.data = diag((-1).^lvec)*gammaLM; % changed sign, now consistent with
% Gamma = gamma0 * exp(i*2*pi*f*(mega*(x2-x1)/c)
glm.lvec = lvec;
glm.mvec = mvec;
glm.Lmax = Lmax;
glm.flow = flow;
glm.deltaF = deltaF;
glm.numFreqs = numFreqs;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a particular gammaLM for debugging
if 1==1
  L=0;
  M=0;

  % find index into gammaLM array  
  indL = find(glm.lvec==L);
  indM = find(glm.mvec==M);
  ii = intersect(indL, indM);

  % construct frequency array
  flow = glm.flow;
  deltaF = glm.deltaF;
  numFreqs = glm.numFreqs;
  f = flow + deltaF*[0:numFreqs-1];
  fhigh = f(end);

  % plot
  plot(f, real(glm.data(ii,:)), f, imag(glm.data(ii,:)));
  axis([flow, fhigh, -.3, .3]);
  grid on;
  xlabel('freq (Hz)');
  ylabel('gammaLM(f)');
  legend(['re, l=',num2str(L),', m=',num2str(M)],...
         ['im, l=',num2str(L),', m=',num2str(M)]);
end

return
