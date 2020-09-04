function [] = sumSums()
% First create sumXFisher files using sumXFisher_condor.
% Then copy the files to dir.
%dir='/archive/home/ethrane/S5H1L1_sumXFisher_L20/';
dir='/archive/home/ethrane/S5H1L1_sumXFisherFlat_L20/';
L=20;
x_tot=zeros((L+1)^2,1);
fisher_tot=zeros((L+1)^2);
coh_tot=0;
segs_tot=0;
for ii=1:18837
  try 
%    load([dir 'sumXFisher_job' num2str(ii) '.mat']);
    load([dir 'sumXFisherFlat_job' num2str(ii) '.mat']);
    x_tot=x_tot+x_opt;
    fisher_tot=fisher_tot+fisher_opt;
    coh_tot=coh_tot+coh;
    segs_tot=segs_tot+segtot;

    % Running Cls-------------------------------
    invfisher_tot=reginv(fisher_tot,12,2/3);
    p_tot=invfisher_tot*x_tot;
    [Cl,ClSig]=constructCl(p_tot,invfisher_tot);
    cl(:,ii)=Cl;
    clSig(:,ii)=ClSig;
  catch
    fprintf('%i\n',ii);
  end
end

%save('/archive/home/ethrane/sumSums.mat');
save('/archive/home/ethrane/sumSumsFlat.mat');
