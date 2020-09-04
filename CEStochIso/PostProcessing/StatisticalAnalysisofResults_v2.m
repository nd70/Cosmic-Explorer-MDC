%function StatisticalAnalysisofResults_v2(slidingOmegasCutCombined, ...
%     concatCCStats, concatNaiveSigmas, dSigmaCut,segmentDuration, ...
%     resampleRate,window1,window2,ifoPair,fileOut,doOverlap,DOFscalefactor, ...
%     outputFileNamePrefix, displayResults, applyBadGPSTimes, printFormats)
function StatisticalAnalysisofResults_v2(slidingOmegasCutCombined, ...
     concatCCStats, concatNaiveSigmas, dSigmaCut,segmentDuration, ...
     resampleRate,window1,window2,ifoPair,fileOut,doOverlap,DOFscalefactor, ...
     outputFileNamePrefix, displayResults, applyBadGPSTimes, printFormats)
% Stochastic background postprocessing script
%
% This script does the following:
%
% 1. Assumes there are two lists of files (concatCCStats and concatNaiveSigmas)
% that correspond to [GPSTimes, CCStatistic, CCsigma] triplets for two types of
% analyses:
%       - naive PSD estimator;
%       - sliding PSD estimator (removes bias of point estimate).
% 
% 2. Using the sigma values for the two types of analysis, the script
% provides graphics to show the presence of outliers in Omega,
% Omega/sigma, etc. due to time-variable PSDs. By suitably editing the
% variable <dSigmaCut>, subset of the full dataset can be analyzed after
% outliers are eliminated.
%
% 3. The final results include a y = A*t + B trend analysis (requires
% curvefit toolbox).
%
% 4. Output data for [time, Omega] values that make the dSigmaCut threshold
% are written to a file for subsequent analysis by the Lombscargle script.
%
% $Id: StatisticalAnalysisofResults_v2.m,v 1.9 2006-09-11 23:18:27 nvf Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyyymmdd = datestr(now,29);

try
  printFormats;
catch
  printFormats = { 'epsc2', 'png' };
end;

% If just a single format, make into a cellarray
if (~iscell(printFormats))
  printFormats = { printFormats };  
end;

% Decide whether to 1. do nothing, 2. calculate & apply delta sigma cuts, or 3.remove user-specified bad GPS times. If 1., then some plots comparing data cuts do not make sense and these will be suppressed.
if applyBadGPSTimes==0
  badTimesCheck=0; %Neither apply user-specified bad GPS times nor apply delta sigma cut on the fly
elseif applyBadGPSTimes==1
  badTimesCheck=1; %Calculate & apply delta sigma cut on the fly
else
  badTimesCheck=2; %User-specified bad GPS times
end

% READ IN THE PIPELINE RESULTS FOR SLIDING PSD ESTIMATES
fileSlidingPSDResults = concatCCStats;
fprintf('Files of sliding results = %s to %s\n', fileSlidingPSDResults{1}, fileSlidingPSDResults{end});
[gpsTimes, stats, sigmas] = readCCStatsFromFile(concatCCStats);
fprintf('Read in %d x 3\n', length(gpsTimes));

% load data into arrays for later manipulation
slidingTimesAll = gpsTimes.';
slidingOmegasAll = stats.'/segmentDuration;
slidingSigmasAll = sigmas.'/segmentDuration;

% READ IN THE PIPELINE RESULTS FOR NAIVE PSD ESTIMATES
if badTimesCheck~=0
  fileNaiveResults = concatNaiveSigmas;
  fprintf('Naive results file(s) = %s ...\n', fileNaiveResults{1});
  [gpsTimes, naiveSigmas, sigmas] = readSigmasFromFile(concatNaiveSigmas);
  fprintf('Read in %d x 3\n',length(gpsTimes));
else
  naiveSigmas = sigmas;
end;

% load data into arrays for later manipulation
naiveTimesAll = gpsTimes.';
naiveSigmasAll = naiveSigmas'./segmentDuration;

% calculate combined quantities for the sliding PSD data
[slidingOmegasAllCombined, slidingSigmasAllCombined, numSegmentsAll] = ...
  combineResultsFromFile('stats', concatCCStats, segmentDuration, ...
                         doOverlap, window1, window2, [], 1);
fprintf('Immediately preceding point estimate, sigma, and number of segments do not include cuts.\n');

% form random deviates for the sliding PSD data
slidingDeviatesAll = (slidingOmegasAll - slidingOmegasAllCombined)./...
                     (slidingSigmasAll);

% find intersection of naive and sliding PSD estimate times for which
% there are both types of values available
[commontimes,slidingIndicesCommon,naiveIndicesCommon] = ...
  intersect(slidingTimesAll,naiveTimesAll);
naiveSigmas = naiveSigmasAll(naiveIndicesCommon);
slidingSigmas = slidingSigmasAll(slidingIndicesCommon);
slidingOmegas = slidingOmegasAll(slidingIndicesCommon);
slidingDeviates = slidingDeviatesAll(slidingIndicesCommon);

% perform outlier cut (trust naive sigmas more)
if badTimesCheck==2
  btimes=load(applyBadGPSTimes);
%  dSigmas = abs(naiveSigmas(ismember(naiveTimesAll,btimes)) - slidingSigmas(slidingTimesAll,btimes))./(naiveSigma(ismember(naiveTimesAll,btimes)));
  dSigmas = abs(naiveSigmas - slidingSigmas)./(naiveSigmas);
  ctimes = find(~ismember(naiveTimesAll,btimes));
  ctimesAll = cumsum(ones(1, length(naiveSigmas)));
else
  dSigmas = abs(naiveSigmas - slidingSigmas)./(naiveSigmas);
  ctimes = find(dSigmas < dSigmaCut);
  ctimesAll = cumsum(ones(1, length(naiveSigmas)));
end

% display three measures representing how much data are rejected by cuts:
fprintf('Effect of various cuts on amount of data analyzed:\n')
if badTimesCheck==2
  fprintf('Number of segments passing user-specified cuts = %d\n', length(ctimes));
  fprintf('Fraction (all sliding PSD data)/(all naive PSD data) = %g\n',...
    length(slidingSigmas)/length(naiveSigmasAll))
  fprintf('Fraction (data passing user-specified cuts)/(all naive PSD data) = %g\n',...
    length(ctimes)/length(naiveSigmasAll))
  fprintf('Fraction (data passing user-specified cuts)/(all sliding PSD data) = %g\n',...
    length(ctimes)/length(slidingSigmasAll))
else
  fprintf('Number of segments passing the dSigmaCut = %d\n', length(ctimes));
  fprintf('Fraction (all sliding PSD data)/(all naive PSD data) = %g\n',...
    length(slidingSigmas)/length(naiveSigmasAll))
  fprintf('Fraction (data passing dSigmaCut)/(all naive PSD data) = %g\n',...
    length(ctimes)/length(naiveSigmasAll))
  fprintf('Fraction (data passing dSigmaCut)/(all sliding PSD data) = %g\n',...
    length(ctimes)/length(slidingSigmasAll))
end

% extract sigmas, omegas, times that satisfy the dSigma cut or user-specified cuts
slidingSigmasCut = slidingSigmas(ctimes);
slidingOmegasCut = slidingOmegas(ctimes);
timesCut = commontimes(ctimes);          
badGPSTimes = transpose(setdiff(slidingTimesAll,timesCut));

% write out bad gps times to file
%save(fileOut, '-ascii','-double','badGPSTimes');
f = fopen(fileOut, 'w');
fprintf(f,'%d\n', badGPSTimes);
fclose(f);

% form random deviates for the sliding PSD data passing the dSigmaCut
slidingDeviatesCut = (slidingOmegasCut - slidingOmegasCutCombined)./...
                     (slidingSigmasCut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
nx = 25; ny = 250;

if (badTimesCheck==1 | badTimesCheck==2)
  set(0,'defaultaxesfontsize',14);
  set(0,'defaulttextfontsize',14);
  h = figure(1);
  %set(h,'Position',[nx,ny,500,500])
  subplot(2,2,1);
  plot(slidingSigmasAll,slidingDeviatesAll,'r.')
  hold on;
  plot(slidingSigmasCut,slidingDeviatesCut,'.')
  xlabel('\sigma');ylabel('(\Omega-<\Omega>)/\sigma');
  title(['Results from ' ifoPair ' for Abs[\delta\sigma]/ \sigma < ' num2str(dSigmaCut)]);
  
  daysSliding = (slidingTimesAll - slidingTimesAll(1))./86400;
  daysCut = (timesCut - timesCut(1))./86400;
  
  subplot(2,2,2);
  plot(daysSliding,slidingOmegasAll,'r');
  hold on;
  plot(daysCut,slidingOmegasCut);
  xlabel('Days since start of run');ylabel('\Omega');
  title(yyyymmdd);
  
  subplot(2,2,3);
  plot(daysSliding,slidingSigmasAll,'r');
  hold on;
  plot(daysCut,slidingSigmasCut);
  xlabel('Days since start of run');ylabel('\sigma');
  hL = legend('All data','Data after Abs[\delta\sigma]/ \sigma cut');
  
  subplot(2,2,4);
  plot(daysSliding,slidingDeviatesAll,'r');
  hold on;
  plot(daysCut,slidingDeviatesCut);
  xlabel('Days since start of run');ylabel('(\Omega-<\Omega>)/\sigma');
  
  filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot1_',ifoPair);
  set(h,'PaperPositionMode','auto','PaperType','usletter',...
    'PaperSize',[8.5 11])
  format compact;
  for pf = printFormats
    print(h, [ '-d' pf{1} ], filename);
  end;
  format;
%  saveas(h,filename,'fig');
  if (~displayResults)
    close(h);
  end;
  % Check whether we're running compiled matlab
  try
    deployedFlag = isdeployed;
  catch
    deployedFlag = 0;
  end;
  if (displayResults && deployedFlag)
    delete(h);
  end;
end; % if (badTimesCheck==1 | badTimesCheck==2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure 2
%h = figure(2); 
%set(h,'Position',[nx+15,ny-15,500,500])
%
%subplot(2,2,1);
%plot(slidingTimesAll,slidingSigmasAll,'r')
%xlabel('GPS Time (s)');ylabel('\sigma');
%set(gca,'Fontsize',14);
%title(['Results from ',ifoPair,' all data']);
%
%subplot(2,2,2);
%plot(slidingTimesAll,slidingDeviatesAll,'r')
%xlabel('GPS Time (s)');ylabel('(\Omega-<\Omega>)/\sigma');
%set(gca,'Fontsize',14);
%title(['Results from ',ifoPair,' all data']);
%
%subplot(2,2,3);
%plot(timesCut,slidingSigmasCut);
%xlabel('GPS Time (s)');ylabel('\sigma');
%set(gca,'Fontsize',14);
%title(['Results from ',ifoPair,' outlier cut']);
%
%subplot(2,2,4);
%plot(timesCut,slidingDeviatesCut);
%xlabel('GPS Time (s)');ylabel('(\Omega-<\Omega>)/\sigma');
%set(gca,'Fontsize',14);
%title(['Results from ',ifoPair,' outlier cut']);
%
%%[max(dSigmas(ctimesAll))
%%[min(dSigmas), max(dSigmas)]
%%[min(naiveSigmas), max(naiveSigmas)]
%
%filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot2_',ifoPair);
%set(h,'PaperPositionMode','auto','PaperType','usletter',...
%  'PaperSize',[8.5 11])
%print(h, '-depsc2',filename)
%print(h, '-dpng',filename)
%saveas(h,filename,'fig');
%if (~displayResults)
%    close(h);
%end;
%if (displayResults && deployedFlag)
%    delete(h);
%end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3
if (badTimesCheck==1 | badTimesCheck==2)
  h = figure(3);
  set(h,'Position',[nx+2*15,ny-2*15,500,500])
  
  subplot(2,1,1);
  [n1,x1] = hist(dSigmas(ctimesAll), [0.0001:1/80.:1]);
  bar(x1,n1,'r');
  hold on;
  [n2,x2] = hist(dSigmas(ctimes), [0.0001:1/80.:1]);
  bar(x2,n2,'c');
  set(gca,'xlim',[0,1]);
  xlabel('Abs[\delta\sigma]/ \sigma');ylabel('# per bin');
  legend('All data','Data after Abs[\delta\sigma]/ \sigma outlier cut');
  title(['Results from ',ifoPair,': Histogram of Abs[\delta\sigma] ' yyyymmdd]);
  
  subplot(2,1,2);
  minx1 = min(slidingSigmas(ctimes)); maxx1 = max(slidingSigmas(ctimes));
  lowx1 = varLimit(minx1,0); highx1 = varLimit(maxx1,1);
  % Hack for some sort of bug...
  if (highx1 < lowx1)
    tmp = lowx1;
    lowx1 = highx1;
    highx1 = tmp;
    clear tmp;
  end;
  nx1 = 50; dx1 = (highx1 - lowx1)/nx1;
  [n1,x1] = hist(slidingSigmas(ctimesAll), [lowx1:dx1:highx1]);
  bar(x1,log10(n1+0.01),'r');
  hold on;
  [n2,x2] = hist(slidingSigmas(ctimes),[lowx1:dx1:highx1]);
  bar(x2,log10(n2+0.01),'c');
  set(gca,'xlim',[0,highx1],'ylim',[0,4]);
  xlabel('\sigma');ylabel('Log10[# per bin]');
  
  filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot3_',ifoPair);
  set(h,'PaperPositionMode','auto','PaperType','usletter',...
    'PaperSize',[8.5 11])
  format compact;
  for pf = printFormats
    print(h, [ '-d' pf{1} ], filename);
  end;
  format;
%  saveas(h,filename,'fig');
  if (~displayResults)
    close(h);
  end;
  if (displayResults && deployedFlag)
    delete(h);
  end;
end; % (badTimesCheck==1 | badTimesCheck==2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4
h = figure(4);
set(h,'Position',[nx+3*15,ny-3*15,500,500])

plot(dSigmas(ctimesAll),naiveSigmas(ctimesAll), 'r.');
hold on;
plot(dSigmas(ctimes), naiveSigmas(ctimes), '.');
xlabel('Abs[\delta\sigma]/\sigma');ylabel('\sigma');
title(['Results from ',ifoPair,': \sigma vs. Abs[\delta\sigma]/ \sigma outlier cut ' yyyymmdd]);

filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot4_',ifoPair);
set(h,'PaperPositionMode','auto','PaperType','usletter',...
  'PaperSize',[8.5 11])
format compact;
for pf = printFormats
  print(h, [ '-d' pf{1} ], filename);
end;
format;
%saveas(h,filename,'fig');
if (~displayResults)
  close(h);
end;
deployedFlag=0; %hack, sGc
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 5 
if (badTimesCheck==1 | badTimesCheck==2)
  h = figure(5);
  set(h,'Position',[nx+4*15,ny-4*15,500,500])
  
  subplot(2,1,1);
  plot(dSigmas, slidingDeviates,'r.');
  hold on;
  plot(dSigmas(ctimes), slidingDeviatesCut,'b.');
  xlabel('Abs[\delta\sigma]/\sigma');ylabel('(\Omega-<\Omega>)/\sigma');
  set(gca,'xlim',[0,2],'ylim',[-20,20]);
  title(['Results from ',ifoPair,': (\Omega-<\Omega>)/\sigma vs. Abs[\delta\sigma]/ \sigma outlier cut ' yyyymmdd]);
  
  subplot(2,1,2);
  plot(slidingSigmas, slidingDeviates,'r.');
  hold on;
  plot(slidingSigmasCut, slidingDeviatesCut,'b.');
  xlabel('\sigma');ylabel('(\Omega-<\Omega>)/\sigma');
  minx1 = min(slidingSigmasCut); maxx1 = max(slidingSigmasCut);
  lowx1 = varLimit(minx1,0); highx1 = varLimit(maxx1,1);
  % Hack for some sort of bug...
  if (highx1 < lowx1)
    tmp = lowx1;
    lowx1 = highx1;
    highx1 = tmp;
    clear tmp;
  end;
  set(gca,'xlim',[0,highx1],'ylim',[-20,20]);
  
  filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot5_',ifoPair);
  set(h,'PaperPositionMode','auto','PaperType','usletter',...
    'PaperSize',[8.5 11])
  format compact;
  for pf = printFormats
    print(h, [ '-d' pf{1} ], filename);
  end;
%  saveas(h,filename,'png');
  format;
%  saveas(h,filename,'fig');
  if (~displayResults)
    close(h);
  end;
  if (displayResults && deployedFlag)
    delete(h);
  end;
end; % if (badTimesCheck==1 | badTimesCheck==2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 6
if (badTimesCheck==1 | badTimesCheck==2)
  h = figure(6);
  set(h,'Position',[nx+4*15,ny-4*15,500,500])
  [n1,x1] = hist(slidingDeviatesAll,[-50:1:50]);
  bar(x1,log10(n1+0.01),'r');
  hold on;
  [n2,x2] = hist(slidingDeviatesCut,[-50:1:50]);
  bar(x2,log10(n2+0.01),'c');
  xlabel('(\Omega-<\Omega>)/ \sigma');ylabel('Log10[# per bin]');
  set(gca,'ylim',[0,floor(max(log10(n1+0.01)))+1]);
  title(['Results from ',ifoPair,': Histogram of (\Omega-<\Omega>)/\sigma with outlier cuts ' yyyymmdd]);
  legend('All data','Data after Abs[\delta\sigma]/ \sigma outlier cut');
  
  filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot6_',ifoPair);
  set(h,'PaperPositionMode','auto','PaperType','usletter',...
    'PaperSize',[8.5 11])
  format compact;
  for pf = printFormats
    print(h, [ '-d' pf{1} ], filename);
  end;
%  saveas(h,filename,'png');
  format;
%  saveas(h,filename,'fig');
  if (~displayResults)
    close(h);
  end;
  if (displayResults && deployedFlag)
    delete(h);
  end;
end; % if (badTimesCheck==1 | badTimesCheck==2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 7
h = figure(7);
set(h,'Position',[nx+5*15,ny-5*15,500,500])

% set up Kolmogorov-Smirnov test on normalized deviates
KSDeviates = sort(slidingDeviatesCut);
ndev  = length(KSDeviates);
KSx = 1/ndev*[1:ndev];% fractional distribution of sorted population (0,...,1)

% Calculate K-S for different sigmas to see dependence and to 
% identify best sigma for data
ddd=[];
sss = [.95:.002:1.1];
for i = 1:length(sss)   
  ddd = [ddd,max(abs(KSx - (erf(((KSDeviates/(sss(i)*sqrt(2))) )) + 1)/2))];
end

plot(sss,ddd,'-r')
xlabel('\sigma')
ylabel('max(abs(data-Erf[...]))')
legend(['Results from ',ifoPair]);
title(['Effect of changing \sigma in Kolmogorov-Smirnov test: max(abs(...)) vs. \sigma ' yyyymmdd],'Fontsize',12)

% use empirically determined best sigma for subsequent analysis
[theMin,indx] = min(ddd);
scaleFactor = sss(indx);

KSPopulationModel = (erf(KSDeviates/(sqrt(2)*scaleFactor))+1)/2;
KSPopulationModelDiff = ...
  (KSx - ( erf(((KSDeviates/(scaleFactor*sqrt(2))) )) + 1)/2);
KSv = max(abs(KSPopulationModelDiff));
nDOF = DOFscalefactor*ndev;
lambda = (sqrt(nDOF) + .12 + .11/sqrt(nDOF))*KSv;
KSValue = KSStatistic(lambda);

fprintf('Best standard deviation describing normalized residuals from KS test = %e\n', scaleFactor);  
fprintf('Minimum of KS test = %e\n', KSv);
fprintf('KS statistic value = %e\n', KSValue);

filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot7_',ifoPair);
set(h,'PaperPositionMode','auto','PaperType','usletter',...
  'PaperSize',[8.5 11])
format compact;
for pf = printFormats
  print(h, [ '-d' pf{1} ], filename);
end;
%saveas(h,filename,'png');
format;
%saveas(h,filename,'fig');
if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 8
h = figure(8);
set(h,'Position',[nx+6*15,ny-6*15,500,500])

subplot(2,1,1);
plot(KSDeviates,KSx,'k.','Markersize',8)
hold on;
plot(KSDeviates,KSPopulationModel,'-c')
axis([-5 5 0 1.1])
title(['Results from ',ifoPair, ': Kolomogorov-Smirnov test for Gaussianity ' yyyymmdd],'Fontsize',12);
legend('Data',['Erf with \sigma=',num2str(scaleFactor,3)],'Location','SouthEast')
text(-4.5,.8,['Kolmogorov-Smirnov test:',num2str(KSValue,2)],'Fontsize',12)
hold off;

subplot(2,1,2)
plot(KSDeviates,KSPopulationModelDiff)
title('Difference between data and Erf(...)','Fontsize',12)
axis([-5 5 -10e-03 10e-03]);
text(-4.5,.005,['Maximum absolute difference:',num2str(KSv,2)],'Fontsize',12)

filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot8_',ifoPair);
set(h,'PaperPositionMode','auto','PaperType','usletter',...
  'PaperSize',[8.5 11])
format compact;
for pf = printFormats
  print(h, [ '-d' pf{1} ], filename);
end;
%saveas(h,filename,'png');
format;
%saveas(h,filename,'fig');
if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 9
h = figure(9);
set(h,'Position',[nx+7*15,ny-7*15,500,500])
dev2 = sum(slidingDeviatesCut.^2);
ndev = length(slidingOmegasCut);
muXY = cov(slidingSigmasCut',slidingDeviatesCut');

spacing = .45;
offset = 0.1;
psize = 0.4;
h1 = subplot('Position',[offset,offset,psize,psize]);
subplot(h1);
minx1 = min(slidingSigmasCut); maxx1 = max(slidingSigmasCut)/4;
lowx1 = varLimit(minx1,0); highx1 = varLimit(maxx1,1);
% Hack for some sort of bug...
if (highx1 < lowx1)
  tmp = lowx1;
  lowx1 = highx1;
  highx1 = tmp;
  clear tmp;
end;

% Occasionally highx1 and lowx can end up being the same
if (lowx1 == highx1)
  highx1 = 10*lowx1;
end;

plot(slidingSigmasCut,slidingDeviatesCut,'.','Markersize',1)
axis([0 highx1 -6 6 ]);
%axis([minx1 max(slidingSigmasCut) min(slidingDeviatesCut) max(slidingDeviatesCut)]);
%keyboard
text(.0,20,['Scatter plot of (\Omega - <\Omega>)/\sigma vs. \sigma for ',ifoPair, ' with outlier cuts ' yyyymmdd],'Fontsize',14)
ylabel('(\Omega - <\Omega>)/\sigma')
xlabel('\sigma')
text(.1*highx1,5,['N = ',num2str(length(slidingSigmasCut)),' pts'],'Fontsize',12)
text(.1*highx1,3.5,['\sigma_{xy}^2 = ',num2str(muXY(1,1))],'Fontsize',12)
text(.1*highx1,2.5,['\rho = \sigma_{xy}^2/\sigma_{x}\sigma_{y}',num2str(muXY(1,2)/sqrt(muXY(1,1)*muXY(2,2)))],'Fontsize',12)
text(.1*highx1,1.5,['\Sigma (\Omega_{i} - <\Omega>)^2/\sigma_{i}^2 = ',num2str(floor(sum(slidingDeviatesCut.^2)))],'Fontsize',12)

h2 = subplot('Position',[offset,spacing+offset/2,psize,psize]);
subplot(h2);
nx1 = 50; dx1 = (highx1 - lowx1)/nx1;
[n1, x1] = hist(slidingSigmasCut, [lowx1:dx1:highx1]);
%stairs(x1, n1 ,'b','Linewidth',1);
stairs(x1, n1);
ym = max(n1);ymax = varLimit(ym,1);
axis([0 highx1 0 ymax ]); %
%axis([minx1 max(slidingSigmasCut) 0 ymax]);
%keyboard
ylabel('Relative Frequency, N')
set(gca,'XTick',[])

h3 = subplot('Position',[spacing+offset/2,offset,psize,psize]);
subplot(h3);
binsize = 0.2;
[n1, x1] = hist(slidingDeviatesCut,[-5:binsize:5]);
muR = mean(slidingDeviatesCut);
sigmaR = scaleFactor;
theory = exp(-(x1 - muR + binsize/2).^2/(2*sigmaR^2));
tmax = 10^log10(max(n1));
ym = max(n1); ymax = varLimit(ym,1);
plot(tmax*theory,x1,'k.','Linewidth',1)
hold on;
%stairs(n1  ,x1 + 0*binsize/2,'b','Linewidth',1)
stairs(n1  ,x1 + 0*binsize/2)
axis([ 0 ymax -6 6]);
%axis([0 ymax min(slidingDeviatesCut) max(slidingDeviatesCut)]);
%keyboard
xlabel('Relative Frequency, N')
set(gca,'YTick',[])
text(.2*ymax,4,['Fit: e^{-x^2/(2\sigma^2)}; \sigma =  ',num2str(sigmaR,'%5.3f')],'Fontsize',12)
text(.2*ymax,3,['Kolmogorov-Smirnov test:',num2str(KSValue,'%4.2f')],'Fontsize',12)

filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot9_',ifoPair);
set(h,'PaperPositionMode','auto','PaperType','usletter',...
  'PaperSize',[8.5 11])
format compact;
for pf = printFormats
  print(h, [ '-d' pf{1} ], filename);
end;
%saveas(h,filename,'png');
format;
%saveas(h,filename,'fig');
if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 11
h = figure(11);
set(h,'Position',[nx+9*15,ny-9*15,500,500])
sigma2Mean = mean(slidingSigmasCut.^2);
hist(1/sigma2Mean*slidingSigmasCut.^2, [0:.05:6])
title(['Results from ', ifoPair, ': Distribution of \sigma^2/<\sigma^2> ' yyyymmdd]);
xlabel('\sigma^2/<\sigma^2>')
ylabel('#/bin')
gpsTimes = timesCut';
relativeTimes = gpsTimes - gpsTimes(1);
filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot11_',ifoPair);
set(h,'PaperPositionMode','auto','PaperType','usletter',...
  'PaperSize',[8.5 11])
format compact;
for pf = printFormats
  print(h, [ '-d' pf{1} ], filename);
end;
%saveas(h,filename,'png');
format;
%saveas(h,filename,'fig');
if (~displayResults)
  close(h);
end;
if (displayResults && deployedFlag)
  delete(h);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uncomment to save valid times and omegas to file for LombScargle analysis
%
%data1 = [relativeTimes,slidingOmegasCut']';
%fileOut = '../matlab/LombScargleAlgorithm/ccstats_FullH1H2_Hann_GPSTimes_Mask5_dSigmaCut_.dat';
%fid = fopen(fileOut,'w')
%fprintf(fid, '%g\t%g\n',data1);
%fclose(fid)

if length(slidingOmegasCut)>3
  days = (relativeTimes)/(3600*24);
  data = slidingOmegasCut';
  errors = slidingSigmasCut';
  opts = fitoptions('Method','LinearLeastSquares','Weights',1./slidingSigmasCut);
  results = fit((days - mean(days))/(days(end) - days(1)),data,'poly1',opts);
  [results, gof,output] = fit((days - mean(days))/(sqrt(var(days))),data,'poly1',opts);
  coeffs = coeffvalues(results);
  confints = confint(results);
  
  h = figure(12);
  set(h,'Position',[nx+10*15,ny-10*15,500,500])
  
  % Plot data
  plot(days,data,'.','Markersize',1)
  hold on
  
  % Plot 3-sigma envelope
  plot(days,3*errors,'b-','Linewidth',2);
  plot(days,-3*errors,'b-','Linewidth',2);
  
  % Plot amplified trend
  plot(days,coeffs(1)*100*((days - (days(end) + days(1))/2)/(days(end) - days(1))) + coeffs(2),'r','Linewidth',1);
  xlabel('Days of Observation')
  ylabel('\Omega_{i}')
  title(['\Omega vs time for each segment ' yyyymmdd])
  scale = sqrt(var(data));
  ymin = -6*scale; ymax = 6*scale;
  lowy1 = varLimit(ymin,0); highy1 = varLimit(ymax,1);
  axis([0 max(days) lowy1 highy1]);
  text(mean(days)/2,0.8*highy1,'Linear trend analysis: \Omega(t) = C_{1}*(t - T_{obs}/2)/T_{obs} + C_{2}','FontSize',10)
  text(mean(days)/2,0.7*highy1,'Line plotted with slope amplified 100X','FontSize',10)
  text(mean(days)/2,0.6*highy1,['C_{1} = ',num2str(coeffs(1),'%6.2e'),'; 95% CB: (',num2str(confints(1,1),'%6.2e'),' , ',num2str(confints(2,1),'%6.2e'),')'],'FontSize',10)
  text(mean(days)/2,0.5*highy1,['C_{2} = ',num2str(coeffs(2),'%6.2e'),'; 95% CB: (',num2str(confints(1,2),'%6.2e'),' , ',num2str(confints(2,2),'%6.2e'),')'],'FontSize',10)
  % write to file
  filename = strcat(outputFileNamePrefix,'StatisticalAnalysis_PanelPlot12_',ifoPair);
  set(h,'PaperPositionMode','auto','PaperType','usletter','PaperSize',[8.5 11])
  for pf = printFormats
    print(h, [ '-d' pf{1} ], filename);
  end;
%  saveas(h,filename,'fig');
  
  if (~displayResults)
    close(h);
  end;
  if (displayResults && deployedFlag)
    delete(h);
  end;
else
  fprintf(1,'Not enough data points to fit a line.  Skipping...\n');
end;
return;
