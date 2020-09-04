% test simulateSB and splice_data routines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some magic values (change as desired)
tlow = 0;
flow = 50;
fhigh = 250;
deltaF = 0.25;
segmentDuration = 60;
resampleRate = 1024;
N = resampleRate*segmentDuration;
deltaT = 1/resampleRate;
numFreqs = floor((fhigh-flow)/deltaF)+1;
signalType = 'const';
detector1 = getdetector('LHO');
detector2 = getdetector('LLO');
H100 = HubbleConstant;
                                                                               
% identity transfer functions
transfer1.flow = flow;
transfer1.deltaF = deltaF;
transfer1.data = ones(numFreqs,1);
transfer2.flow = flow;
transfer2.deltaF = deltaF;
transfer2.data = ones(numFreqs,1);
                                                                                
% initialize arrays for time-series data
tseries1 = [];
tseries2 = [];

for k=1:10

  % simulate SB signals
  [h1_left, h2_left]   = simulateSB(tlow, deltaT, deltaT, N, N, ...
                                    signalType, ...
                                    detector1, detector2, ...
                                    transfer1, transfer2, ...
                                    0, 0, 0, 0, NaN, NaN);
                                                                                
  [h1_right, h2_right] = simulateSB(tlow, deltaT, deltaT, N, N, ...
                                    signalType, ...
                                    detector1, detector2, ...
                                    transfer1, transfer2, ...
                                    0, 0, 0, 0, NaN, NaN);
                       
  % splice the data                                                         
  data1 = spliceData(1, transpose(h1_left.data), transpose(h1_right.data), ...
                     N, 1);
  data2 = spliceData(1, transpose(h2_left.data), transpose(h2_right.data), ...
                     N, 2);

  % concatenate with previous portion of the signal
  tseries1 = [tseries1 data1];
  tseries2 = [tseries2 data2];

  % construct time series data
  h1 = constructTimeSeries(tseries1, tlow, deltaT);
  h2 = constructTimeSeries(tseries2, tlow, deltaT);

  % discrete times (assumed to be consistent between h1 and h2)
  t = h1.tlow + h1.deltaT*[0:length(h1.data)-1]';                                                                               
  % plot time series
  figure(1);
  subplot(2,1,1); plot(t, h1.data);
  title('Time series','FontSize',12)
  ylabel('Amplitude 1','FontSize',12);
                                                                                
  subplot(2,1,2); plot(t, h2.data);
  xlabel('time (sec)','FontSize',12);
  ylabel('Amplitude 2','FontSize',12);
  
  % calculate power spectra              
  [temp1, f] = psd(real(h1.data),[],1/deltaT);
  [temp2, f] = psd(real(h2.data),[],1/deltaT);
  P1 = 2*deltaT*temp1;
  P2 = 2*deltaT*temp2;
                                                                                
  % plot power spectra
  figure(2);
  subplot(2,1,1); loglog(f, P1);
  title('Power spectra','FontSize',12);
  ylabel('Log10(P1) (strain^2/Hz)','FontSize',12);
                                                                                
  subplot(2,1,2); loglog(f, P2);
  xlabel('Log10(frequency) (Hz)','FontSize',12);
  ylabel('Log10(P2) (strain^2/Hz)','FontSize',12);
                            
end
                   
return
