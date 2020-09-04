function [ out ] = plotPositionSlice( position, t, Xf, f, mask, res ,F)
% plots a time shift slice at some position
%
% inputs:
%  position: [Nx2] array, 1st column: RA, 2nd c.: DECL
%  t:   time vector
%  Xf:  map output [(Lmax+1)^2,fmin...fmax,-fmin...-fmax]
%  f:   frequency mask
%  mask: freqency mask
% 
%

global PIXELCONVERSION;

N=size(Xf,1);
ff=[f(:);-f(:)];
Lmax=sqrt(N)-1;
try
 res; 
 catch 
 res=1; 
 end
checkPixelConversion(Lmax,res);

Usub=zeros(size(position,1),N);
for ii=1:N
    Usub(:,ii)=interpolateMap( PIXELCONVERSION.U(:,ii), ...
                                     PIXELCONVERSION.RA, ...
                                     PIXELCONVERSION.DECL, position );
end

try
 
    VV=real(Usub*F*Usub');
    ss=sqrt(diag(VV))*ones(1,length(t));

 catch 

    ss=ones(size(position,1),length(t));
end

xx=real(Usub*Xf*exp(i*2*pi*ff*t)) ./ss;
plot(t,xx);
grid on;
xlabel('time shift [sec]')
ylabel('X')
for jj=1:size(position,1)
    leg{jj}=sprintf('RA%.1f DECL%.1f',position(jj,1),position(jj,2));
end
legend(leg);

if nargout>0
    out=xx;
end

end

