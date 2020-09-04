function [invfisher,pCovar,egmod,U,eg,V,pCovar_UL]=reginv(fisher, regMethod, regCutoff)
% invert fisher matrix after regularization
% Parameters:
%   input:
%    fisher:     Fisher matrix to be inverted
%    regMethod:  one of the following values
%                 1: Keep all eigenvalues with eg(ii)/max(eg) >= regCutoff
%                 2: Keep a fraction of all eigenvalues, regCutoff is that fraction; 
%                 3: Keep N=regCutoff of all eigenvalues; 
%                    1-3: All modified eigenvalues are set to the smallest non-modified value. 
%                11: Keep all eigenvalues with eg(ii)/max(eg) >= regCutoff
%                12: Keep a fraction of all eigenvalues, regCutoff is that fraction; 
%                13: Keep N=regCutoff of all eigenvalues; 
%                    11-13: All modified eigenvalues are set to infinity, i.e. corresponding
%                           eigenvalue of invfisher is set to zero.  
%    regCutoff: cut off value, see above
%               Default is regMethod=2; regCutoff=2/3;
%
%   output:
%    invfisher:  Inverted Fisher matrix (after regularisation)
%    pCovar:     covariance of the regularized result, pCovar = V*diag(eg./egmod.^2)*U';
%    egmod:      modified eigenvalues of the Fisher matrix, i.e. invfisher=V*diag(1./egmod)*U'
%    eg:         original eigenvalues of Fisher matrix, i.e. fisher=U*diag(eg)*V'
%    U, V:       unitary matrices from SVD
%    pCovar_UL:  modified covariance of the regularized result,
%                used for UL calculation. It assumes that the average
%                power in the invisible modes is equal to the power in the
%                used modes:
%                  eg_UL = eg;
%                  eg_UL(ii+1:end)=1/mean(1./eg(1:ii));
%                where ii is the last kept eigenmode, and
%                  pCovar_UL = V*diag(1./eg_UL)*U';
%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stefan Ballmer, sballmer@caltech.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  regMethod;
catch
  regMethod = 2;
end;
try
  regCutoff;
catch
  regCutoff = 2/3;
end;

[U,S,V] =svd(fisher);
eg=diag(S);
egmod=eg;

switch regMethod
  case {1,2,3}
    mult=1;
  case {11,12,13}
    mult=0;
  otherwise,
    warning('Unknown regularizarion method');
end;    

switch regMethod
  case {1,11} 
    ind=find( eg./eg(1) >= regCutoff);
    if length(ind)>0, ii=ind(end); else ii=1; end;
  case {2,12}
    ii=floor(length(eg)*regCutoff);
  case {3,13}
    ii=regCutoff;
  otherwise,
    warning('Unknown regularizarion method.');
end;    
if(ii<1), ii=1; end;
if(ii>length(eg)), ii=length(eg); end;

if mult==1
  egmod(ii+1:end)=eg(ii);
else
  egmod(ii+1:end)=Inf;
end
invfisher = V*diag( 1./egmod   )*U';
pCovar    = V*diag(eg./egmod.^2)*U';

% if required, calculate the extended covariance for the upper limit
% calculation
if nargout==7
    eg_UL = eg;
    eg_UL(ii+1:end)=1/mean(1./eg(1:ii));
    pCovar_UL = V*diag(1./eg_UL)*U';
end

