function [varargout]=randomTS(dur,fsample,Hf,intLog,Npol,tau,NFFT,varargin)
%  function [h1,h2,...]=randomTS(dur,fsample,Hf,intLog,Npol,tau,NFFT,transfer1,transfer2,...)
%
%  randomTS  -- creates nargout independent random time series 
%               with one-sided power spectrum Hf/Npol
%                    
%  arguments: dur     - duration of the time series in seconds
%             fsample - sample frequency of the time series
%             Hf      - total power spectrum (one-sided) in all polarizations, i.e.
%                       the actual power spectrum of each time series 
%                       hi is Hf / Npol
%                       preferred input: Nx2 array (filename and freq. series
%                       work too ... in principle)
%             intLog  - boolean whether to interpolate Hf logarithmicly
%             Npol    - number of polarizations
%             tau     - vector of time delays, produces hi of the form
%                       [hi(t-tau1),hi(t-tau2),...]
%                       argument optional, default 0
%             NFFT    - N for the fft (optional, smallest possible is used if ommited)
%             transfer- transfer functions to calibrate detector i
%
%  output:    hi   - data of the ith time series
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some fixed parameters:
% number fo points for the hanning  window smoothening of the calibration cutt-off
Nhann=1000;
% Nhann corresponds to dFhann Hz
dFhann=20;

Nsample=ceil(fsample*dur);

try
  tau;
catch
  tau = 0;
end;
try
  NFFT;
catch
  NFFT = 2^ceil(log2(Nsample));
end;

if length(NFFT) == 0
  NFFT = 2^ceil(log2(Nsample));
end;

if and(nargin > 7,nargin ~= 7+length(tau))
  error('Number of time delays tau have to agree with number of transfer functions');
end

% make sure NFFT is big enough
if NFFT < Nsample
  error('NFFT < Nsample');   
end

%numFreqs=length(Hf.data);
%flow=Hf.flow;
%deltaF=Hf.deltaF;
numFreqs=NFFT/2;
flow=0;
deltaF=fsample/NFFT;
f=(0:1:(numFreqs-1))'*deltaF+flow;

nlow=flow/deltaF;
nhigh=nlow+numFreqs-1;

fdoublesided=fftshift([0,(-NFFT/2+1):(NFFT/2-1)])'*deltaF;
for kk=1:length(tau)
  if tau(kk)~=0
    shfter(:,kk)=exp(i*2*pi*tau(kk).*fdoublesided);
  end;
end;

H=checkArgumentFreqSeries(Hf,flow,deltaF,numFreqs,intLog);
% normalization: H/2/Npol to get 2-sided PSD, 1/sqrt(2) for complex, deltaF for the integration
nnorm=sqrt(H.data.*(deltaF/(2*2*Npol)))*NFFT;

% interpolate transfer functions to appropriate frequenciex
if nargin > 7
  hn=hann(Nhann);
  for kk=1:length(tau)
    transferarg=varargin{kk};
    d_transfer=transferarg.data;
    f_transfer   = transferarg.flow + transferarg.deltaF*[0:length(transferarg.data)-1]';
    % add hanning portions for smooth cutt-off
    if transferarg.flow>0
      % end of hann portion fits the transfer function
      hnup=hn(1:(Nhann/2))*d_transfer(1);
      % if enough room below flow: use full dFhann/2, otherwise use whatever is left
      if transferarg.flow > dFhann/2
        fhnup=((-Nhann/2):(-1))'*(dFhann/Nhann) + transferarg.flow;
      else
        fhnup=(0:(Nhann/2-1))'*(2*transferarg.flow/Nhann);
      end;
      if fhnup(1)<=0
        f_transfer=[fhnup;f_transfer];
        d_transfer=[hnup; d_transfer];
      else
        f_transfer=[0;fhnup;f_transfer];
        d_transfer=[0;hnup; d_transfer];
      end;
    end
    fhigh=f(end);
    tfhigh=f_transfer(end);
    if tfhigh < fhigh
      % end of hann portion fits the transfer function
      hndn=hn((1+Nhann/2):(Nhann))*d_transfer(end);
      % if enough room above f_transfer(end): use full dFhann/2, otherwise use whatever is left
      if fhigh-tfhigh > dFhann/2
        fhndn=(1:(Nhann/2))'*(dFhann/Nhann) + tfhigh;
      else
        fhndn=(1:(Nhann/2))'*(2*(fhigh-tfhigh)/Nhann) + tfhigh;
      end;
      if fhndn(Nhann/2)>=fhigh
        f_transfer=[f_transfer;fhndn];
        d_transfer=[d_transfer;hndn];
      else
        f_transfer=[f_transfer;fhndn;fhigh];
        d_transfer=[d_transfer;hndn;0];
      end;
    end
    transfer{kk} = interp1(f_transfer, d_transfer, f);
  end;
end;

% loop through all hi
for ii=1:nargout
  nbase = nnorm .* ( randn(numFreqs, 1) + i* randn(numFreqs, 1) );
  for kk=1:length(tau)
    if nargin > 7
      %calibrate
      n=nbase.*transfer{kk};
    else
      n=nbase;
    end;
    nflip=fliplr(flipud(conj(n)));
    % start filling in hf, leave DC and Nyquest = 0
    hf= zeros(NFFT,1);
    if nlow>0
      hf([nlow:nhigh]+1)=n;
      hf([(NFFT-nhigh):(NFFT-nlow)]+1)=nflip;
    else
      hf([nlow+1:nhigh]+1)=n(2:numFreqs);
      hf([(NFFT-nhigh):(NFFT-nlow-1)]+1)=nflip(1:(numFreqs-1));
    end;
    if tau(kk)==0
      h(:,kk)=real(ifft(hf));
    else
      h(:,kk)=real(ifft(hf.*shfter(:,kk)));
    end;
  end;
  %crop output 
  varargout(ii)={h(1:Nsample,:)};
end;
return;
