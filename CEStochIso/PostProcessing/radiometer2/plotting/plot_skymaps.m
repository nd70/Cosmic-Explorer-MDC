function plot_skymaps(pproc_params)

if isstr(pproc_params)
       load(pproc_params);
end
%load('O1Week1Radiometer_20_297-skymap_COMBINED-JOBS-0-1-163');
TAG = getFileTag(pproc_params, [1:pproc_params.numJobs]);
load([pproc_params.directory '/' pproc_params.prefix '_COMBINED-JOBS' TAG]);
SNR = FINAL_COMBINED.pte.data ./ FINAL_COMBINED.sigma.data;

skyPattern = SkyPattern(360,181);
r=unique(skyPattern(:,1));
d=unique(skyPattern(:,2));
K=length(r);
L=length(d);
ra=reshape(skyPattern(:,1),L,K);
decl=reshape(skyPattern(:,2),L,K);
picSNR = reshape(SNR,L,K);
picSigma = reshape(FINAL_COMBINED.sigma.data, L, K);
picY = reshape(FINAL_COMBINED.pte.data, L, K);


figure;
plot_map(picSNR, ra*360/24, decl, 'SNR Map','SNR');
print('-dpng','-r1000',[pproc_params.directory '/' pproc_params.prefix 'SNR_MAP' TAG]);
print('-depsc2',[pproc_params.directory '/' pproc_params.prefix 'SNR_MAP' TAG]);
close;

figure;
plot_map(log10(picSigma), ra*360/24, decl, 'SIGMA','SNR');
print('-dpng','-r1000',[pproc_params.directory '/' pproc_params.prefix 'SIGMA_MAP' TAG]);
close;

figure;
plot_map(log10(abs(picY)), ra*360/24, decl, 'POINT ESTIMATE','SNR');
print('-dpng','-r1000',[pproc_params.directory '/' pproc_params.prefix 'POINT_ESTIMATE_MAP' TAG]);
print('-depsc2',[pproc_params.directory '/' pproc_params.prefix 'POINT_ESTIMATE_MAP' TAG]);
close;

function plot_map(pic, ra, dec, plot_title, cbarlabel)
m_proj('hammer','lon',[0 360],'lat',[-90 90]);
m_pcolor(ra, dec, pic);
shading flat;
m_grid;
xlabel('\phi');
ylabel('\theta');
h = title(plot_title);
P = get(h,'Position');
set(h,'Position',[P(1) P(2)+0.3 P(3)]);
% Set colorbar and its title, and adjust positions.
c = colorbar('location','SouthOutside');
P1 = get(c,'Position');
set(c,'Position',[P1(1) P1(2)-0.13 P1(3) P1(4)]);
set(get(c, 'Xlabel'), 'String', cbarlabel);
% Adjust color axis scale.
%caxis([-1 1]);

