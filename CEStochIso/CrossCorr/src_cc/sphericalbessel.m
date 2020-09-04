function y = sphericalbessel(n,z)
%  SPHERICALBESSEL -- Calculate spherical Bessel functions
%
%  sphericalbessel(n,z) returns values of spherical Bessel function 
%  j_n(z), order n evaluated at z
%
%  Any zero elements of z have the corresponding output elements set
%  to z.^n (Matlab interprets 0^0 as 1); the results for non-zero 
%  arguments are calculated by calling BESSELJ with a fractional
%  order and multiplying by the appropriate factor, i.e., sqrt(pi/(2z))
% 
%  $Id: sphericalbessel.m,v 1.1 2007-04-10 19:34:47 jromano Exp $
%  Routine adapted by John T. Whelan from one written by Joe Romano.
%  Contact john.whelan@ligo.org
% 
%  See also BESSELJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default value for zero elements of z
y = z.^n;

% correct value for non-zero elements of z
i = find(z);
y(i) = ((0.5*pi./z(i)).^0.5).*besselj(n+.5, z(i)); % see page 437 Abramowitz

return;

