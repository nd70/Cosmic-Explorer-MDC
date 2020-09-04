function gamma = sidtimeavgdcorf(det1, det2, decDeg);
%  SIDTIMEAVGDCORF.M -- Calculate DC point sorce overlap reduction function
%                              averaged over siderial time
%
%  sidtimeavgdcorf.m(f, det1, det2, decDeg) calculates the overlap
%  reduction function for gravitational waves coming from a pointlike
%  unpolarized source at the specified declination decDeg (given in
%  degrees), averaged over siderial time, for the detectors det1 and
%  det2.  Both detectors are structures in the standard Cartesian
%  detector geometry format.  This structure has the fields
%     r: [3x1 double] % position vector (in units of meters)
%                       in Earth-based Cartesian coordinates
%     d: [3x3 double] % response tensor in Earth-based Cartesian coordinates
% 
%  Given geographical location information, detector geometry structures
%  can be built using the functions BUILDIFODETECTOR and
%  BUILDBARDETECTOR.
% 
%  Note that the calculation of this code is done independently of that
%  in ORFINTEGRAND, but analytically speaking, the former result
%  should be obtainable from the latter by averaging over right
%  ascension or hour angle.
% 
%  Note that since this is the average of gamma, rather than the
%  root-mean-square, it's not actually all that useful.
% 
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
% 
%  See also BUILDIFODETECTOR, BUILDBARDETECTOR, CREATEDETECTOR,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of inputs

error(nargchk(3,3,nargin));

LAL_C_SI = 2.99792458e8; % m/s

% trace-free response tensors
d1 = det1.d - eye(3).*trace(det1.d)/3;
d2 = det2.d - eye(3).*trace(det2.d)/3;

% calculate c1, c2, c3 coeffs.
c1 = sum(sum(d1.*d2));
c2 = d1(3,:) * d2(:,3);
c3 = d1(3,3) * d2(3,3);

cosdelta2 = (cos(decDeg*pi/180)).^2 ;
cosdelta4 = cosdelta2.^2;

zeta1 =   1.0  -       cosdelta2 + 0.125  * cosdelta4;
zeta2 = - 2.0  + 4.0 * cosdelta2 - 1.25   * cosdelta4;
zeta3 =   0.5  - 2.5 * cosdelta2 + 2.1875 * cosdelta4;

gamma = 5 * (c1*zeta1 + c2*zeta2 + c3*zeta3);

return;
