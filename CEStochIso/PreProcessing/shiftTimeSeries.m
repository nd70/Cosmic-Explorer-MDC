function h_shift=shiftTimeSeries(h0,fsample,tau)
%shifts time series of hplus and hcross to take into account the time shift
% tau between the two detectors.

NFFT=length(h0);
h0_tilde=fft(h0)/NFFT;
deltaF=fsample/NFFT;
temp=fftshift([0,(-NFFT/2+1):(NFFT/2-1)])'; 
fdoublesided=temp*deltaF;

%debugging info for identifying interesting frequencies for narrow band injections
%max_h0_tilde=max(abs(h0_tilde));
%relevant_cut=(abs(h0_tilde)>1e-4*max_h0_tilde);


%multiply h0_tilde by exp(2pi*i*f*tau) to transform h0(t) to h0(t+tau)
Tf=exp(2*pi*i*fdoublesided*tau);
h0_tilde_shift=Tf.*h0_tilde;


%ifft the new h0_tilde to get the new h
h_shift=ifft(h0_tilde_shift)*NFFT;


return

