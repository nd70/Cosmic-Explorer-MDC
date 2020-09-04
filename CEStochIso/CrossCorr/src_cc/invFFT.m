function [xt,varargout]=invFFT(xf,tmax,NFFT)
%  function [xt,NFFT]=invFFT(xf,tmax,NFFT)
%
%  Calculates xt(t) = ind_-Inf^Inf df xf(f) exp(i*2*pi*f*t)
%
%  input:  xf:   frequency series (not heterodyned, symmetry=1) 
%                also provides flow, deltaF, and numFreqs
%          tmax: maximum that is needed
%          NFFT: N for the FFT (optional, default = 32 times oversampled)
%
%  output: xt:   time series,
%                deltaT=1/(deltaF * NFFT)
%                length(xt)=floor(deltaF * NFFT * tmax)
%                The data is filled in starting at tlow (which is negative),
%                i.e. the zero shift entry is xt(length(xt)/2+1).
%          NFFT  returns the used NFFT (optional)
%
%                the output is calculated via
%
%                xt(n) = 2 x deltaF x
%                       Re[ exp(i*2*pi*flow*n/(deltaF*NFFT) x 
%                           NFFT x ifft(xf) ] 
%
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaulfOversample=32;
numFreqs=length(xf.data);
flow=xf.flow;
deltaF=xf.deltaF;
fhigh=flow+deltaF*(numFreqs-1);
numhigh=fhigh/deltaF;
if xf.symmetry == 0
  error('no heterodyned data support in invFFT');
end;
try
  NFFT;
catch
  NFFT = defaulfOversample * 2^ceil(log2(numhigh));
end;
if length(NFFT)==0
  NFFT = defaulfOversample * 2^ceil(log2(numhigh));
end
if NFFT < numFreqs
  error('NFFT too small');
end;
nmax=floor(deltaF * NFFT * tmax);
fsample=deltaF*NFFT;
deltaT=1/fsample;

n=(-nmax:nmax-1)';
t=n*deltaT;
tlow=t(1);

norm=2*deltaF*NFFT;
xftilde=fftshift(ifft(xf.data,NFFT));
if flow == 0
  xt.data=norm*real(xftilde(1+NFFT/2+n))-deltaF*real(xf.data(1));
else
  xt.data=norm*real(exp(i*2*pi*flow*t).*xftilde(1+NFFT/2+n));
end;
xt.tlow=tlow;
xt.deltaT=deltaT;
xt.fbase=NaN;
xt.phase=NaN;
if nargout==2
  varargout(1)={NFFT};
end
return;

