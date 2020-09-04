function standard_pproc_plots(pproc_params, job1, njobs, badjobs)

if isstr(pproc_params)
    load(pproc_params);
end
try 
    badjobs;
catch
    badjobs = [];
end

LB=flipud(lbmap(256,'BrownBlue'));

LOADTAG = getFileTag(pproc_params, [1:pproc_params.numJobs]);
TAG = getFileTag(pproc_params,[job1:job1+njobs-1]);
jobs_to_plot = [job1:job1+njobs-1];
load([pproc_params.directory '/' pproc_params.prefix '_INDIVIDUAL-JOBS' LOADTAG]);
load([pproc_params.directory '/' pproc_params.prefix '_COMBINED-JOBS' LOADTAG]);


RESIDUALS = zeros(length(jobs_to_plot), length(SINGLEJOBS{1}.pte.f));
SIGMAS = zeros(length(jobs_to_plot), length(SINGLEJOBS{1}.pte.f));

bins = [-5:0.1:5];
RESIDUALS_HIST = zeros(length(bins), length(SINGLEJOBS{1}.pte.f))';

for ii = 1:length(jobs_to_plot)
    i = jobs_to_plot(ii);
    if any(badjobs == i)
        RESIDUALS(ii,:) = NaN*ones(1, length(SINGLEJOBS{1}.pte.f));
        SIGMAS(ii,:) = NaN*ones(1, length(SINGLEJOBS{1}.pte.f));
        RESIDUALS_HIST(ii,:) = NaN*ones(1,length(bins));
        fprintf('HI THERE!!!')
        continue
    end
    try
    RESIDUALS(ii,:) = (SINGLEJOBS{i}.pte.data - FINAL_COMBINED.pte.data) ./ SINGLEJOBS{i}.sigma.data;
    SIGMAS(ii,:) = SINGLEJOBS{i}.sigma.data;
    catch
    RESIDUALS(ii,:) = NaN*ones(1, length(SINGLEJOBS{1}.pte.f));
    SIGMAS(ii,:) = NaN*ones(1, length(SINGLEJOBS{1}.pte.f));
    RESIDUALS_HIST(ii,:) = NaN*ones(1,length(bins));
    end
end
ADTESTS = zeros(length(RESIDUALS_HIST(1,:)), 1);
for i = 1:length(SINGLEJOBS{1}.pte.f)
    idxs = find( RESIDUALS(:,i) == 0 ); 
    RESIDUALS(idxs,1) = NaN;
    RESIDUALS_HIST(i,:) = histc(RESIDUALS(~isnan(RESIDUALS(:,i)),i), bins);
    if or(isempty(find(RESIDUALS(:,i))), SINGLEJOBS{1}.pte.data(i)==0)
        ADTESTS(i) = NaN;
        continue
    end
    try
    [h1,ADTESTS(i)] = adtest(RESIDUALS(~isnan(RESIDUALS(:,i)),i));
    catch
    ADTESTS(i) = NaN;
    end
end
figure;
printmap(RESIDUALS', jobs_to_plot, SINGLEJOBS{1}.pte.f, 'Job Number','Frequency [Hz]', 'RESIDUAL',[-5,5]);
colormap(LB), h=colorbar;
title('map of residuals');
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_INDIVIDUAL-JOBS-RESIDUALS' TAG])
close;
SIGMAS(SIGMAS==inf) = NaN;


figure;
printmap(log10(SIGMAS)', jobs_to_plot, SINGLEJOBS{1}.pte.f, 'Job Number','Frequency [Hz]', 'SIGMA');
colormap(parula);
h=colorbar;
title('map of sigmas');
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_INDIVIDUAL-JOBS-SIGMAS' TAG])
close;

figure;
printmap(log10(RESIDUALS_HIST), bins, SINGLEJOBS{1}.pte.f, 'Job Number','Frequency [Hz]', 'RESIDUAL',[0,log10(15)]);
colormap(LB), h2=colorbar;
title('map of residuals');
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_INDIVIDUAL-JOBS-RESIDUALS-HISTOGRAM' TAG])
close;
figure;
scatter(SINGLEJOBS{1}.pte.f,ADTESTS);
title('map of residuals');
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_ADTESTS-SCATTER' TAG])
close;

figure;
hist(ADTESTS)
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_ADTESTS-PDF' TAG])
close;
