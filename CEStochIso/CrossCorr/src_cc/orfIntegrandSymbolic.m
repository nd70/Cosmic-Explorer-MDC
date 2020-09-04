function [gamma]=orfIntegrandSymbolic(det1,det2,time,ra,decl)
%  function [gamma]=orfIntegrandSymbolic(det1,det2,time,ra,decl)
%
%  orfIntegrandSymbolic  -- calculates the projection F(ra,decl,pol)
%                           on to detector det at sidereal time time
%                           Output is a struct with gamma0, tau and
%                           F1p,F1x,F2p,F2x
%                    
%  orfIntegrand:
%             det1 - detector 1 structure containing position r and tensor d
%             det2 - detector 2 structure containing position r and tensor d
%             time - sidereal time in hours (0..24h)
%             ra   - right ascension in hours of source in the sky
%             decl - right ascension in degrees of source in the sky
%
%  averages over polarizations
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% separation vector between sites
s = det2.r - det1.r;

% distance between sites
distance = norm(s);

% unit vector Omega
ez=[0;0;1];
Omega=rotateVector(ez,pi/12*(ra-time),pi/2-pi/180*decl,0);

% time delay
c=299792458;
tau = dot(Omega,s)/c;

gamma.F1p=antenna(det1,time,ra,decl,   0);
gamma.F1x=antenna(det1,time,ra,decl,pi/4);
gamma.F2p=antenna(det2,time,ra,decl,   0);
gamma.F2x=antenna(det2,time,ra,decl,pi/4);

% divide by 2 since Hf is defined as sum over pol of power
gamma.gamma0=(gamma.F1p * gamma.F2p + ...
             gamma.F1x * gamma.F2x)/2;
	     
gamma.tau=tau;
