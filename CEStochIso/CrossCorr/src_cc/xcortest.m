  resampleRate1=4096;
  resampleRate2=4096;
  fsample=resampleRate1;
  flow=70;
  deltaF=0.25;
  numFreqs=7721;
  f=(0:(numFreqs-1))'*deltaF+flow;
  alphaBetaFile1='/home/sballmer/stochastic/sgwb/S4/input/calibration/UnityAlphaBetas.mat';
  calCavGainFile1='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-H1-DERRCAL-CAV_GAIN-793099715.mat';
  calResponseFile1='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-H1-CAL-DERRRESPONSE-793099715.mat';
  ASQchannel1='LSC-DARM_ERR';
  alphaBetaFile2='/home/sballmer/stochastic/sgwb/S4/input/calibration/UnityAlphaBetas.mat';
  calCavGainFile2='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-L1-DERRCAL-CAV_GAIN-792579613.mat';
  calResponseFile2='/home/sballmer/stochastic/sgwb/S4/input/calibration/S04-L1-CAL-DERRRESPONSE-792579613.mat';
  %calResponseFile2=calResponseFile1; calCavGainFile2=calCavGainFile1;
  ASQchannel2='LSC-DARM_ERR';
  [t1, f1, R01, C01, alpha1, gamma1] = ...
    readCalibrationFromFiles(alphaBetaFile1, calCavGainFile1, calResponseFile1);
  [t2, f2, R02, C02, alpha2, gamma2] = ...
    readCalibrationFromFiles(alphaBetaFile2, calCavGainFile2, calResponseFile2);
  
  H1=load(calResponseFile1);
  L1=load(calResponseFile2);
  H1.R0(1)=Inf;
  L1.R0(1)=Inf;
  
  %skyPattern=load('/home/sballmer/stochastic/sgwb/S4/matlab/apps/FullSkyPattern_360_x_181.txt');
  skyPattern=SkyPattern(360,181);

  detector1=getdetector('LHO');
  detector2=getdetector('LLO');
  filename='/home/sballmer/stochastic/sgwb/S4/input/radiometer/sim/pointSources.txt';
  simPSSkyPositions=load(filename);
  filename='/home/sballmer/stochastic/sgwb/S4/input/radiometer/sim/Hf.txt';
  simPSPowerSpec=load(filename);
  Hf=checkArgumentFreqSeries(simPSPowerSpec,flow,deltaF,numFreqs);
  % prepare the point source simulation
  segmentDuration=60;
  bufforget=512; % just has to be long enough
  bufspline=32; % power of 2 makes it more efficient
  ra=simPSSkyPositions(:,1);
  decl=simPSSkyPositions(:,2);
  if size(simPSSkyPositions,2)>2
    power=simPSSkyPositions(:,3);
  else
    power=ones(size(ra));
  end
  initPointSourceData(resampleRate1,simPSPowerSpec,true,detector1,detector2,ra,decl,power,flow,deltaF,numFreqs,...
                      t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
		      t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
		      bufforget,bufspline);


pspec=1;

bufferSecs1=1;
bufferSecs2=1;

hannDuration1 = 60;
numPoints1    = segmentDuration*resampleRate1; 
dataWindow1   = tukeywin(numPoints1, hannDuration1/segmentDuration);
fftLength1    = 2*numPoints1;

hannDuration2 = 60;
numPoints2    = segmentDuration*resampleRate2;
dataWindow2   = tukeywin(numPoints2, hannDuration2/segmentDuration);
fftLength2    = 2*numPoints2;



% begin
GPS=793197092;
dur=segmentDuration+2*bufferSecs1;
GPSmid=GPS+dur/2;
sid=GPStoGreenwichMeanSiderealTime(GPS);
sidmid=GPStoGreenwichMeanSiderealTime(GPSmid);
g=orfIntegrandSymbolic(detector1,detector2,sidmid,ra,decl);

% calib
        calibsec1 = GPS + bufferSecs1;
        [R1, responseOK1] = ...
          calculateResponse(t1, f1, R01, C01, alpha1, gamma1, calibsec1,...
			    ASQchannel1);

        % if response function is bad, set flag and exit loop
        if responseOK1==false
          badResponse = true;
          break
        end
  
        % evaluate response function at desired frequencies
        response1 = convertResponse(f1, R1, flow, deltaF, numFreqs, 0, 0);

        % convert to transfer functions (units: counts/strain) 
        transfer1 = convertResponse(f1, R1, flow, deltaF, numFreqs, 1, 0);
        calibsec2 = GPS + bufferSecs2;
        [R2, responseOK2] = ...
          calculateResponse(t2, f2, R02, C02, alpha2, gamma2, calibsec2,...
			    ASQchannel2);

        % if response function is bad, set flag and exit loop
        if responseOK2==false
          badResponse = true;
          break
        end
   
        % evaluate response function at desired frequencies
        response2 = convertResponse(f2, R2, flow, deltaF, numFreqs, 0, 0);

        % convert to transfer functions (units: counts/strain) 
        transfer2 = convertResponse(f2, R2, flow, deltaF, numFreqs, 1, 0);
  


[h1, h2] = getPointSourceData(GPS,dur);
o1 = constructTimeSeries(h1, ...
                         GPS,  1/fsample, ...
			 NaN, 0);
o2 = constructTimeSeries(h2, ...
                         GPS,  1/fsample, ...
			 NaN, 0);


if false
  Ndata=11;
  for q=1:(Ndata-1)
    G=GPS+q*dur;
    [h1a, h2a] = getPointSourceData(G,dur);
    o1 = constructTimeSeries(o1.data+h1a, ...
                             GPS,  1/fsample, ...
			     NaN, 0);
    o2 = constructTimeSeries(o2.data+h2a, ...
                             GPS,  1/fsample, ...
			     NaN, 0);  
  end
  o1 = constructTimeSeries(o1.data/Ndata, ...
                           GPS,  1/fsample, ...
			   NaN, 0);
  o2 = constructTimeSeries(o2.data/Ndata, ...
                           GPS,  1/fsample, ...
			   NaN, 0);    
end

        highpassed1 = o1;
        highpassed2 = o2;

    p=ones(numFreqs,1);
    calPSD1_avg = constructFreqSeries(p, flow, deltaF, 1);
    calPSD2_avg = constructFreqSeries(p, flow, deltaF, 1);

    mask=constructFreqSeries(p, flow, deltaF, 1);

      firstIndex1 = 1 + bufferSecs1*resampleRate1;
      firstIndex2 = 1 + bufferSecs2*resampleRate2;
      lastIndex1  = length(highpassed1.data)-bufferSecs1*resampleRate1;
      lastIndex2  = length(highpassed2.data)-bufferSecs2*resampleRate2;

      r1 = constructTimeSeries(highpassed1.data(firstIndex1:lastIndex1), ...
                                  highpassed1.tlow + bufferSecs1, ...
                                  highpassed1.deltaT, ...
				  highpassed1.fbase, highpassed1.phase);
      r2 = constructTimeSeries(highpassed2.data(firstIndex2:lastIndex2), ...
                                  highpassed2.tlow + bufferSecs2, ...
                                  highpassed2.deltaT, ...
				  highpassed2.fbase, highpassed2.phase);

    rbartilde1 = windowAndFFT(r1, dataWindow1, fftLength1);
    rbartilde2 = windowAndFFT(r2, dataWindow2, fftLength2);

    gamma=constructFreqSeries(ones(numFreqs,1),flow,deltaF,1);
    [Q, ccVar, sensInt] = calOptimalFilter(segmentDuration, gamma, ...
                                           Hf, [], ... 
                                           calPSD1_avg, calPSD2_avg, ...
                                           dataWindow1, dataWindow2, mask);

      [ccStat,ccSpec] = calCrossCorr(rbartilde1, rbartilde2, Q, ...
                                   response1, ...
                                   response2,0.011);
				   
      fprintf('Starting readout\n');
      tic;
      map=ccStatReadout(detector1, detector2, GPSmid, skyPattern, ccStat, ccVar);
      toc;
      			   
      tauccStat=-((0:length(ccStat.data)-1)*ccStat.deltaT + ccStat.tlow);



close all;

tau=min(tauccStat):1e-6:max(tauccStat);
taushift=g.tau;
taux=tau-taushift;
ind=find(taux==0);
if length(ind) >0
  taux(ind)=1;
end;
th=(sin(2*pi*f(end)*taux)-sin(2*pi*f(1)*taux))./(pi*taux);
if length(ind) >0
  th(ind)=2*(f(end)-f(1));
end
th=th/2/1930*60;
tau=taux+taushift;
plot(tau,th,tauccStat,ccStat.data/g.gamma0,tauccStat,ones(size(tauccStat))*sqrt(ccVar),tauccStat,-ones(size(tauccStat))*sqrt(ccVar));
grid;

if 1
  r=unique(skyPattern(:,1));
  d=unique(skyPattern(:,2));
  K=length(r);
  L=length(d);
  pic=zeros(K,L);
  kk=0;
  SNR=map(:,1)./map(:,2);
  fprintf('Start rearranging data\n');
  tic;
  ra=reshape(skyPattern(:,1),L,K);
  decl=reshape(skyPattern(:,2),L,K);
  pic=reshape(SNR,L,K);
  toc;
  
  pic(104,188)=120;
  %figure;
  %ugly but impossible to screw up axis
  %pcolor(ra,decl,pic)
  
  % nice picture - but the axis!
  figure;
  imagesc(r,d,flipud(pic));
  colorbar;
  ft=flipud(get(gca,'YTickLabel'));
  set(gca,'YTickLabel',ft);
  xlabel('Right ascension [hours]');
  ylabel('Declination [deg]');
  title('Radiometer antenna pattern');

end






if pspec
[p1,f1]=pwelch(h1/sqrt((g.F1p^2+g.F1x^2)/2),[],[],[],fsample);
[p2,f2]=pwelch(h2/sqrt((g.F2p^2+g.F2x^2)/2),[],[],[],fsample);
cal1=interp1(H1.f,abs(H1.R0),f1);
cal2=interp1(L1.f,abs(L1.R0),f2);
%cal1=interp1(ftrans,abs(trans1.data),f1);
%cal2=interp1(ftrans,abs(trans2.data),f2);
%cal1=1; cal2=1;
%loglog(f1,p1.*cal1.^2,f2,p2.*cal2.^2,H1.f,abs(H1.R0).^2,L1.f,abs(L1.R0).^2);
%loglog(f1,pp1.*cal1.^2,f1,p1.*cal1.^2,f2,p2.*cal2.^2);
ind1=find(and(f1>100,f1<1000));
ind2=find(and(f2>100,f2<1000));
cp1=p1.*cal1.^2;
cp2=p2.*cal2.^2;
%loglog(f1(ind1),cp1(ind1),f2(ind2),cp2(ind2));
figure;
loglog(f1,cp1,f2,cp2);
grid;
end
%save out3.mat
