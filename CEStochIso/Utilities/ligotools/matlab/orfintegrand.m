function gamma = orfintegrand(f, det1, det2, n);
%  ORFINTEGRAND -- Calculate overlap reduction function integrand
%                              between two detectors
%
%  orfintegrand(f, det1, det2, n) calculates the values of the angular
%  integrand of the overlap reduction function for gravitational waves
%  coming from the direction specified by the unit vector n at the
%  frequencies included in the vector f for the detectors det1 and
%  det2.  Both detectors are structures in the standard Cartesian
%  detector geometry format.  This structure has the fields
%     r: [3x1 double] % position vector (in units of meters)
%                       in Earth-based Cartesian coordinates
%     d: [3x3 double] % response tensor in Earth-based Cartesian coordinates
%  
%  The unit vector n is in the same Earth-fixed Cartesian direction,
%  which will depend on siderial time for a source specified by right
%  ascension and declination.
% 
%  Given geographical location information, detector geometry structures
%  can be built using the functions BUILDIFODETECTOR and
%  BUILDBARDETECTOR.
% 
%  Note that the calculation of this code is done independently of that
%  in OVERLAPREDUCTIONFUNCTION, but analytically speaking, integrating
%  the former result should give the latter.
% 
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
% 
%  See also BUILDIFODETECTOR, BUILDBARDETECTOR, CREATEDETECTOR,
%  GETCARTESIANDIRECTIONFROMSOURCE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of inputs

error(nargchk(4,4,nargin));

% make sure freq vector the right way round
f = f(:);

LAL_C_SI = 2.99792458e8; % m/s

% trace-free response tensors
d1 = det1.d - eye(3).*trace(det1.d)/3;
d2 = det2.d - eye(3).*trace(det2.d)/3;

% Matrices to promote to frequency and direction arrays
profreq = ones(size(f));
sizedir = size(n);
prodir = ones(1,sizedir(2));

% calculate c1, c2, c3 coeffs.
d1n = d1 * n;
d2n = d2 * n;

c1 = sum(sum(d1.*d2)) * prodir;
c2 = sum(d1n .* d2n);
c3 = sum(n .* d1n) .* sum(n .* d2n);

magnitude = profreq * ( 5/(4*pi) * ( c1 - 2 * c2 + c3 / 2 ) );
phase = (2 * pi / LAL_C_SI) * f * ( transpose(det2.r-det1.r) * n ) ;

gamma = magnitude .* exp(1i * phase);

return;
