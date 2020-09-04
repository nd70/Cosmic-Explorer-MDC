function [Cl_UL, C] = cl_limits(sphmat, Ntrials)
% function cl_limits(sphmat, Ntrials)
% calculate upper limits on angular power distribution weights Cl
% 
% sphmat  : mat file containing sph post processing results
% Ntrials : number of trials for upper limits calculation
%
% Cl_UL   : Cl upper limit array
% C       : array containing distribution of each Cl, calculated numerically for Ntrial independant trials
%
% Author  : E Thrane, L. Sammut 
% Updated : Sep 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = Ntrials;                   % number of trials for upper limits

% prepare invfisher as in constructCl_real.m (p0 will serve as plm)
invfisher = sphmat.invfisher0*sphmat.fisher0*sphmat.invfisher0;

% change to real spherical harmonic basis
[C2R, R2C] = getC2RMatrix(sphmat.L);
invfisherR = C2R*invfisher*R2C;

% diagonalize Fisher matrix and get eigenvalues
[V,D] = eig(invfisherR);
egnvls = diag(D);

% get distributions of Cls numerically
for ii=1:T
  z0 = randn(length(invfisherR),1);
  zz = sqrt(egnvls).*z0;     % scale random numbers
  % transform back to sph basis but...
  % ...no need to worry about real part because we use abs(p).^2
  p = V*zz;

  % calculate distribution of Cls
  for l=0:sphmat.L
    ind = getLind_real(sphmat.L, l);
    C(ii,l+1) = sum(abs(p(ind)).^2) - trace(abs(invfisherR(ind,ind)));
    C(ii,l+1) = C(ii,l+1)/(2*l+1);
  end
end

% calcluate upper limit
Cl_UL = getUpperLimitFromMC(C, sphmat.Cl, sphmat.Cl_sig, sphmat.conf, sphmat.calErr);

if any(Cl_UL==0)
  fprintf('error\n');
end

end
