function Mout=invLMblock(M)

s=size(M);

if s(1)~=s(2)
  disp('not a square matrix!');
  return;
end

lmax=sqrt(s(1))-1;

Mout=zeros(size(M));
ii=0;
for m=-lmax:lmax
    ind=ii+(1:lmax-abs(m)+1);
    ii=ii+(lmax-abs(m)+1);
    Mout(ind,ind)=inv(M(ind,ind));
end

