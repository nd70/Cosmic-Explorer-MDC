function hh=plotMapAitoff(map,Nra,Ndecl,outputFilePrefix,spt,colorflg,circra,circdecl)
try spt; catch spt=0; end;
try colorflg; catch colorflg=0; end;
try circra;
  try circdecl; catch
      try
          circdecl=circra.decl; circra=circra.ra;
      catch
          circdecl=circra(:,2);
          circra=circra(:,1);
      end
  end;
end;


hh=0;

  nRows=size(map,1);
  nPoints=Nra*Ndecl;
  % check dimension
  % if just one line too big: that's the isotropic result
  % just ignore it.
  if nRows == nPoints+1
    map=map(1:nPoints,:);
  elseif nRows ~= nPoints
    error('map has wrong dimensions');
    return;
  end
  skyPattern=SkyPattern(Nra,Ndecl);
  r=unique(skyPattern(:,1));
  d=unique(skyPattern(:,2));
  K=length(r);
  L=length(d);
  ra=reshape(skyPattern(:,1),L,K);
  decl=reshape(skyPattern(:,2),L,K);

  sptall=spt;
for spt=sptall
  if spt>=0
    S=map(:,1);
    N=map(:,2);
    SNR=S./N;
    mm=max(SNR);
    mmi=min(SNR);
    vv=sqrt(var(SNR,1));
    fprintf('SNR standard deviation: %g\tMax SNR: %g\tMin SNR: %g\n',vv,mm,mmi);

    picS=reshape(S,L,K);
    picN=reshape(N,L,K);
    pic =reshape(SNR,L,K);
  else
    pic =reshape(map(:,-spt),L,K);
  end;
if spt==1 || spt==0
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  hh=m_pcolor(ra/24*360,decl,fliplr(pic));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18');%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
  bmin=min(-1,min(SNR));
  bmax=max(+1,max(SNR));
  caxis([bmin,bmax]);
  %title('SNR');
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  title(cbr,'SNR');%,'fontangle','italic');
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'snr_map']);
end
if spt==2 || spt==0
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  hh=m_pcolor(ra/24*360,decl,fliplr(picS));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18);%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
  bmin=min(0,min(min(picS)));
  bmax=max(0,max(max(picS)));
  caxis([bmin,bmax]);
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  title('point estimate');%,'fontangle','italic');
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'pte_map']);
end
if spt==3 || spt==0
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  hh=m_pcolor(ra/24*360,decl,fliplr(picN));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18);%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
  bmin=(min(min(picN)));
  bmax=(max(max(picN)));
  caxis([bmin,bmax]);
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  title('theoretical sigma');%,'fontangle','italic');
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'sig_map']);
end
if spt==4 || spt==4
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  picmll=-loglikelihood(pic);
  hh=m_pcolor(ra/24*360,decl,fliplr(picmll));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18);%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
  bmin=0; %(min(min(picmll)));
  bmax=(max(max(picmll)));
  caxis([bmin,bmax]);
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  title('-log(p(x>abs(snr)))');%,'fontangle','italic');
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'll_map']);

end
if spt>=50
  if spt>=5000
    % dirty quick trick: spt=calError+confidence*10000
    confidence=floor(spt)/10000;
    calError=spt-10000*confidence;
    Y   =reshape(picS,L*K,1);
    YSig=reshape(picN,L*K,1);
    syscorr=load('/home/sballmer/stochastic/sgwb/S4/matlab/apps/ttsystematicErrorIncrease_current.mat');
    YSig=YSig.*syscorr.ttRelError(1:nPoints);
    fprintf('Included systematic error\n');
    mu_UL=getUpperLimit(Y,YSig,confidence,calError);
    mu_UL=reshape(mu_UL,L,K);
    fprintf('Confidence Level      : %g%%\n',confidence*100);
    fprintf('Calibration Error     : %g%%\n',calError*100);
  elseif spt>=500
    % dirty quick trick: spt=calError+confidence*1000
    confidence=floor(spt)/1000;
    calError=spt-1000*confidence;
    Y   =reshape(picS,L*K,1);
    YSig=reshape(picN,L*K,1);
    mu_UL=getUpperLimit(Y,YSig,confidence,calError);
    mu_UL=reshape(mu_UL,L,K);
    fprintf('Confidence Level      : %g%%\n',confidence*100);
    fprintf('Calibration Error     : %g%%\n',calError*100);
  else  
    confidence=spt/100; % spt is now confidence level in %
    [mode, mu_UL, bayesfactor] =bayesianUL(picS,picN,confidence);
    fprintf('Confidence Level      : %g%%\n',confidence*100);
  end
  maxUL=max(max(mu_UL));
  minUL=min(min(mu_UL));
  maxInd=find(maxUL==mu_UL);
  minInd=find(minUL==mu_UL);
  fprintf('Maximal upper limit   : %g @ ra= %f, decl= %f\n',maxUL,ra(maxInd),decl(maxInd));
  fprintf('Minimal upper limit   : %g @ ra= %f, decl= %f\n',minUL,ra(minInd),decl(minInd));
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  hh=m_pcolor(ra/24*360,decl,fliplr(mu_UL));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18');%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
  bmin=min(0,min(min(mu_UL)));
  bmax=max(0,max(max(mu_UL)));
  caxis([bmin,bmax]);
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  title([num2str(confidence*100),'% confidence upper limit']);%,'fontangle','italic');
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'upper_limit_map']);
end

if spt<0
  smartfigure;
  if colorflg
    oldmap=colormap;
    bw=(1-(0:63)/63)'*[1,1,1];
    colormap(bw);
  end;
  m_proj('hammer','lon',[0,360],'lat',[-90,90]);
  hh=m_pcolor(ra/24*360,decl,fliplr(pic));
  set(hh,'Linestyle','none');
    m_gridAstro(...
           'xtick',[0,90,180,270,360],...
	   'xticklabels',[24,18,12,6,0],...
	   'ytick',[-90,-45,0,45,90],...
	   'yticklabels',[-90,-45,0,45,90],...
	   'XaxisLocation','middle',...
	   'YaxisLocation','left',...
	   'color',[0,0,0],...
	   'fontsize',25,...
	   'linewidth',4,...
	   'linestyle','-',...
	   'fontname','Helvetica',...
	   'box','off');
xlabel('Right ascension [hours]','fontsize',18);
y=ylabel('Declination [degree]','fontsize',18);%,'fontangle','italic');
  pos=get(y,'position');
  pos(1)=pos(1).*1.05;
  set(y,'position',pos);
% EHT - hack this to force bmin to be zero.
  bmin=(min(min(pic)));
%  bmin=0;
  bmax=(max(max(pic)));
  caxis([bmin,bmax]);
  cbr=colorbar('EastOutside');
%  set(cbr,'DataAspectRatio',[1,0.2*(bmax-bmin),1]);
  try circra; m_range_ring(360-circra/12*180,circdecl,1000); end;
  set(gcf,'Renderer','painters');
  print('-dpng',[outputFilePrefix 'snr_map']);
end

end
try circra;
  [pte,sig,snr]=readOneDirection(map,Nra,Ndecl,circra,circdecl);
  fprintf('Position:        \t\tRA= %g\tDECL= %g\n',circra,circdecl);
  fprintf('Point Esitmate / Sigma: \t%g  +- %g\n',pte,sig);
  fprintf('Signal-to-Noise ratio: \t\t%g\n',snr);
end

set(gca,'FontSize',18); %added by Eric Thrane

return
