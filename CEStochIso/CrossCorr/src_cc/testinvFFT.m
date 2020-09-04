close all;
clear;

f=(00:0.25:2000)';
source='/home/sballmer/stochastic/sgwb/S4/input/radiometer/sourceDirection.txt';
HS='/home/sballmer/stochastic/sgwb/S4/input/radiometer/Hf.txt';
target=[0 45];
target='/home/sballmer/stochastic/sgwb/S4/input/radiometer/sourceDirection.txt';

H ='/home/sballmer/stochastic/sgwb/S4/input/radiometer/Hf.txt';
Hg='/home/sballmer/stochastic/sgwb/S4/input/radiometer/Hgaussian.txt';


calH1='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-H1-CAL-DERRRESPONSE-793099715.mat';
calL1='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-L1-CAL-DERRRESPONSE-792579613.mat';

s=load(source);



Hf=load(H);

H1=load(calH1);
L1=load(calL1);
H1.R0(1)=Inf;
L1.R0(1)=Inf;

flow=f(1);
deltaF=0.25;
fhigh=f(end);
numFreqs=length(f);

dur=2;
fsample=4096;
det1=getdetector('LHO');
det2=getdetector('LLO');
ra=[s(1,1)];
decl=[s(1,2)];
time=0;

[hp,hx]=randomTS(dur,fsample,Hf,true,2,[0,0]);

xf=constructFreqSeries(hp(1:numFreqs,1)+i*hx(1:numFreqs,1),flow,deltaF,1);
tmax=0.01;


[xt,N]=invFFT(xf,tmax,2^15);
t=0:xt.deltaT:(tmax-xt.deltaT);
t=[(flipud(fliplr(-t))-xt.deltaT),t];
fd=unique([-f;f]);
NN=length(fd);

clear xfd
fourierkernel=exp(i*2*pi*fd*t);
xfd(1:numFreqs,1)=flipud(fliplr(conj(xf.data)));
xfd(NN-numFreqs+1:NN,1)=(xf.data);
xtc=deltaF*real(transpose(fourierkernel)*xfd);
norm(xt.data-xtc)/sqrt(norm(xt.data)*norm(xtc))
