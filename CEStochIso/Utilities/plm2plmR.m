function [ plmR ] = plm2plmR( plm )
%PLM2PLMR Summary of this function goes here
% This function does the same as
%   C2R  = getC2RMatrix(lmax);
%   plmR = real(C2R*plm);
% But it uses the old function plm2plmreal, and then resorts the entries
% from the old real Spherical Harmonics format to the new real Spherical
% Harmonics format.
%
% See getC2RMatrix.m for a description of the new and old format.
%
% Written by Stefan Ballmer
% sballmer@ligo.caltech.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lmax=sqrt(length(plm))-1;
    plmreal=plm2plmreal(plm);
    
    index=plmreal(:,1:2);
    data =transpose(plmreal(:,3:4));
    ind  =index~=0;
    ind(1,1)= true;
    indT=transpose(ind);
    
    plmR=reshape(data(indT),(lmax+1)^2,1);
    
end

