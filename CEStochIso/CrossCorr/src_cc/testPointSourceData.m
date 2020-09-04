
sources=load('/home/sballmer/stochastic/sgwb/S4/input/radiometer/sim/pointSources.txt');

ASQchannel1     ='LSC-DARM_ERR';
alphaBetaFile1  ='/home/sballmer/stochastic/sgwb/S4/input/calibration/UnityAlphaBetas.mat';
calCavGainFile1 ='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-H1-DERRCAL-CAV_GAIN-793099715.mat';
calResponseFile1='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-H1-CAL-DERRRESPONSE-793099715.mat';
ASQchannel2     ='LSC-DARM_ERR';
alphaBetaFile2  ='/home/sballmer/stochastic/sgwb/S4/input/calibration/UnityAlphaBetas.mat';
calCavGainFile2 ='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-L1-DERRCAL-CAV_GAIN-792579613.mat';
calResponseFile2='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-L1-CAL-DERRRESPONSE-792579613.mat';

fsample=4096;
Hf=load('/home/sballmer/stochastic/sgwb/S4/input/radiometer/sim/Hf.txt');
det1=getdetector('LHO');
det2=getdetector('LLO');
ra=sources(:,1);
decl=sources(:,2);
  if size(sources,2)>2
    power=sources(:,3);
  else
    power=ones(size(ra));
  end
flow=65;
deltaF=0.25;
numFreqs=7741;
bufforget=256;
bufspline=16;
[t1, f1, R01, C01, alpha1, gamma1] = ...
  readCalibrationFromFiles(alphaBetaFile1, calCavGainFile1, calResponseFile1);
[t2, f2, R02, C02, alpha2, gamma2] = ...
  readCalibrationFromFiles(alphaBetaFile2, calCavGainFile2, calResponseFile2);

% initialize the global variable POINT_SOURCE_MEMORY
initPointSourceData(fsample,Hf,true,det1,det2,ra,decl,power,flow,deltaF,numFreqs,...
                             t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
			     t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
			     bufforget,bufspline);

% everything is ready - start testing

GPSstart=795010000;
dur=30;
dg=dur;

[h1, h2] = getPointSourceData(GPSstart,dg);


Nt=10;
% look-back test
if true
 t=GPSstart;
 [h1i, h2i] = getPointSourceData(t,dg);
 ll=0;
 while ll<Nt
  fprintf('test run %d\n',ll);
  p= round(rand*(256));
  t=t+p;
  [h1n, h2n] = getPointSourceData(t,dg);
  if p<=bufforget-dur-bufspline
    [h1m, h2m] = getPointSourceData(t-p,dg);
    if norm(h1m-h1i)+norm(h2m-h2i) > 0
      error('data not identical');
    else
      fprintf('data identical (look-back: %d sec)\n',p);
      ll=ll+1;
    end
  end;
  h1i=h1n;
  h2i=h2n;
 end
end



N=length(h1(:,1));
% x-correlation test
if false
 Han=h1(1:N).*hann(N)*2;
 Liv=h2(1:N).*hann(N)*2;
 hanf=fft(Han);
 livf=fft(Liv);
 df=fsample/N;
 ind=fftshift([0,(-N/2+1):(N/2-1)]);
 fd=(ind.*df)';
 tau=[-0.013:0.0001:0.013];
 sh=exp(-i*2*pi*fd*tau);
 hl=conj(hanf).*livf;
 r=transpose(hl)*sh;

 plot(tau,abs(r));
 grid;
 figure;
end;

%power spectrum test
if true
 [p1,f1]=pwelch(h1,[],[],[],fsample);
 [p2,f2]=pwelch(h2,[],[],[],fsample);

 calH1='../../input/calibration/S04-H1-CAL-DERRRESPONSE-793099715.mat';
 calL1='../../input/calibration/S04-L1-CAL-DERRRESPONSE-792579613.mat';

 H1=load(calH1);
 L1=load(calL1);
 H1.R0(1)=Inf;
 L1.R0(1)=Inf;
 cal1=interp1(H1.f,abs(H1.R0),f1);
 cal2=interp1(L1.f,abs(L1.R0),f2);
 %cal1=1; cal2=1;
 time=GPStoGreenwichMeanSiderealTime(GPSstart+dur);
 g=orfIntegrandSymbolic(det1,det2,time,ra,decl);
 %[(g.F1p^2+g.F1x^2)/(g.F2p^2+g.F2x^2) , time,ra(1),decl(1)]
 ccc1=2/abs(g.F1x*g.F1x+g.F1p*g.F1p);
 ccc2=2/abs(g.F2x*g.F2x+g.F2p*g.F2p);
 %ccc1=1; ccc2=1;
 loglog(f1,p1.*cal1.^2*ccc1,f2,p2.*cal2.^2*ccc2);
 grid;
end
