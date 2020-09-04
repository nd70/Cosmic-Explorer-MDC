function [invGamma]=invertFisher(Gamma)
%
% inverts fisher information matrix 

% first complete the matrix
%
sz=size(Gamma);
szglm=min(sz);
Lmax=round((2*szglm+.25)^0.5 - 1.5);

mvec=zeros(szglm,1);
lvec=zeros(szglm,1);
for mm=0:Lmax
  i0=mm*(Lmax+1)-mm*(mm-1)/2 +1;
  mvec(i0:(i0+Lmax-mm))=mm;
  lvec(i0:(i0+Lmax-mm))=mm:Lmax;
end
lvecfull=[lvec; lvec(2+Lmax:end)];
mvecfull=[mvec;-mvec(2+Lmax:end)];




if sz(1) > sz(2)
  lmlm=(mvecfull+lvecfull)*ones(1,length(lvec)) + ones(length(lvecfull),1)*(transpose(mvec+lvec));
  Gamma_m=conj(Gamma).*(1-mod(lmlm,2)*2);
  Gamma_m=[Gamma_m(1:1+Lmax,:);Gamma_m(sz(2)+1:end,:);Gamma_m(2+Lmax:sz(2),:)];
  Gamma=[Gamma,Gamma_m(:,2+Lmax:end)];
elseif sz(1) < sz(2)
  lmlm=(mvec+lvec)*ones(1,length(lvecfull)) + ones(length(lvec),1)*(transpose(mvecfull+lvecfull));
  Gamma_m=conj(Gamma).*(1-mod(lmlm,2)*2);
  Gamma_m=[Gamma_m(:,1:1+Lmax),Gamma_m(:,sz(2)+1:end),Gamma_m(:,2+Lmax:sz(2))];
  Gamma=[Gamma;Gamma_m(2+Lmax:end,:)];
end
invGamma=inv(Gamma);
