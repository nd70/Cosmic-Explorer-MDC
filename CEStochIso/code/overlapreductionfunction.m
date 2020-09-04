function gamma = overlapreductionfunction(f, det1, det2);
%  OVERLAPREDUCTIONFUNCTION -- Calculate overlap reduction function 
%                              between two detectors
%
%  overlapreductionfunction(f, det1, det2) calculates the values of the overlap
%  reduction function at the frequencies included in the vector f for
%  the detectors det1 and det2.  Both detectors are structures in the
%  standard Cartesian detector geometry format.  This structure has the
%  fields
%     r: [3x1 double] % position vector (in units of meters)
%                       in Earth-based Cartesian coordinates
%     d: [3x3 double] % response tensor in Earth-based Cartesian coordinates
% 
%  Given geographical location information, detector geometry structures
%  can be built using the functions BUILDIFODETECTOR and
%  BUILDBARDETECTOR.
% 
%  Routine adapted by John T. Whelan from one written by Joe Romano.
%  Contact john.whelan@ligo.org
% 
%  See also BUILDIFODETECTOR, BUILDBARDETECTOR, CREATEDETECTOR,
%  SPHERICALBESSEL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of inputs

error(nargchk(3,3,nargin));

LAL_C_SI = 2.99792458e8; % m/s

% separation vector between sites

s = det2.r - det1.r;

% distance between sites
distance = norm(s);

% trace-free response tensors
d1 = det1.d - eye(3).*trace(det1.d)/3;
d2 = det2.d - eye(3).*trace(det2.d)/3;

if (distance == 0)
  gamma = ones(size(f)) .* 2 * trace(d1*d2);
  return;
end

% unit separation vector
s = s./distance;

% calculate c1, c2, c3 coeffs.
c1 = sum(sum(d1.*d2));
c2 = (transpose(s)) * d1 * d2 * s;
c3 = ((transpose(s)) * d1 * s)*((transpose(s)) * d2 * s);

alpha = f.*(distance*2*pi/LAL_C_SI);

% calculate overlap reduction function
alpha2 = alpha.*alpha;

% Default values for zero arguments
b0 = ones(size(alpha));
b1 = b0 / 3;
b2 = b0 / 15;

% Pick out the non-zero elements of alpha

i = find(alpha);

b0(i) = sphericalbessel(0,alpha(i));
b1(i) = sphericalbessel(1,alpha(i))./alpha(i);
b2(i) = sphericalbessel(2,alpha(i))./alpha2(i);

rho1 =   5.0*b0 - 10.0*b1 +  5.0*b2;
rho2 = -10.0*b0 + 40.0*b1 - 50.0*b2;
rho3 =   2.5*b0 - 25.0*b1 + 87.5*b2;

gamma = c1*rho1 + c2*rho2 + c3*rho3;

return;
