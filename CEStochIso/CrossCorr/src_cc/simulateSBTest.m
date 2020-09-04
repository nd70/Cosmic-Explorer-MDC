function simulateSBTest
%
%  simulateSBTest --- test the simulateSB routine for some simple cases
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
% 
%  $Id: simulateSBTest.m,v 1.3 2005-02-24 16:12:05 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some magic values (change as desired)
seed = 314159; 
tlow = 0;
flow = 50;
fhigh = 250;
deltaF = 0.25;
segmentDuration = 60;
resampleRate = 1024;
N=resampleRate*segmentDuration;
deltaT = 1/resampleRate;
numFreqs = floor((fhigh-flow)/deltaF)+1;
%signalType = 'white';
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

% simulate signal
[h1, h2] = simulateSB(tlow, deltaT, deltaT, N, N, ...
                      signalType, ...
                      detector1, detector2, ...
                      transfer1, transfer2, ...
                      0, 0, 0, 0, NaN, NaN, seed);

% check tlow
fprintf('tlow of h1 time-series = %e\n', h1.tlow);
fprintf('tlow of h2 time-series = %e\n', h2.tlow);

% check deltaT
fprintf('deltaT of h1 time-series = %e\n', h1.deltaT);
fprintf('deltaT of h2 time-series = %e\n', h2.deltaT);

% discrete times (assumed to be consistent between h1 and h2)
t = h1.tlow + h1.deltaT*[0:length(h1.data)-1]';

% check variance
fprintf('Variance of h1 time-series = %e\n', var(real(h1.data)));
fprintf('Variance of h2 time-series = %e\n', var(real(h2.data)));

% plot time series
figure(1);
plot(t,real(h1.data));
title('Time series','FontSize',14);
xlabel('time (sec)','FontSize',14);
ylabel('amplitude','FontSize',14);

figure(2);
plot(t,real(h2.data));
title('Time series','FontSize',14);
xlabel('time (sec)','FontSize',14);
ylabel('amplitude','FontSize',14);

% plot power spectra
[temp1, f] = psd(real(h1.data),[],1/deltaT);
[temp2, f] = psd(real(h2.data),[],1/deltaT);
P1 = 2*deltaT*temp1;
P2 = 2*deltaT*temp2;

figure(3);
plot(f, log10(P1));
title('Power spectrum','FontSize',14);
xlabel('frequency (Hz)','FontSize',14);
ylabel('Log10(Power) (strain^2/Hz)','FontSize',14);

figure(4);
plot(f, log10(P2));
title('Power spectrum','FontSize',14);
xlabel('frequency (Hz)','FontSize',14);
ylabel('Log10(Power) (strain^2/Hz)','FontSize',14);

figure(5);
loglog(f, P1);
title('Power spectrum','FontSize',14);
xlabel('Log10(frequency) (Hz)','FontSize',14);
ylabel('Log10(Power) (strain^2/Hz)','FontSize',14);

figure(6);
loglog(f, P2);
title('Power spectrum','FontSize',14);
xlabel('Log10(frequency) (Hz)','FontSize',14);
ylabel('Log10(Power) (strain^2/Hz)','FontSize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare value of power spectra with expected value at f0

f0=100;

for k = 1:10

  % simulate signal for different seeds k
  [h1, h2] = simulateSB(tlow, deltaT, deltaT, N, N, ...
                        signalType, ...
                        detector1, detector2, ...
                        transfer1, transfer2, ...
                        0, 0, 0, 0, NaN, NaN, seed);


  % calculate power spectra
  [temp1, f] = psd(real(h1.data),[],1/deltaT);
  [temp2, f] = psd(real(h2.data),[],1/deltaT);
  P1 = 2*deltaT*temp1;
  P2 = 2*deltaT*temp2;

  ind=find(f==f0);
  P1_f0(k) = P1(ind);
  P2_f0(k) = P2(ind);

  % display value
  fprintf('trial = %d, P1 at f_0=%d Hz: %e strain^2/Hz\n',k,f0,P1_f0(k));
  fprintf('trial = %d, P2 at f_0=%d Hz: %e strain^2/Hz\n',k,f0,P2_f0(k));
 
end

% display simulated and expected values
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
fprintf('Mean value of P1 at f_0=%d Hz: %e strain^2/Hz\n', f0, mean(P1_f0));
fprintf('Mean value of P2 at f_0=%d Hz: %e strain^2/Hz\n', f0, mean(P2_f0));
fprintf('Stddev of P1 values =%e strain^2/Hz\n', std(P1_f0));
fprintf('Stddev of P2 values =%e strain^2/Hz\n', std(P2_f0));

expected = (3*H100^2)/(10*pi^2) * (1/f0^3);
fprintf('Expected value of power at f_0=%d Hz: %e strain^2/Hz\n', f0, expected);

return

