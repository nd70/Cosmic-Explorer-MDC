function plot_psd_map_single_job(pproc_params,ifonum,jobnum)

if isstr(pproc_params)
    load(pproc_params);
end
[psdjob] = read_psd_file(pproc_params.paramsFile, ifonum, jobnum)
DATA = zeros(size(psdjob.data));
for i = 1:length(psdjob.f)   
    DATA(i,:) = psdjob.data(i,:);
    idxs = (DATA(i,:) / median(psdjob.data(i,:)))> 7;
    DATA(i,idxs) = NaN;
end

printmap(log10(psdjob.data), psdjob.times, psdjob.f, 'Time','Frequency [Hz]', ['PSD ' num2str(ifonum)]);
colormap(parula), h=colorbar;
title(['Map of PSD ' num2str(ifonum)]);
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_INIDVIDUAL-JOB-' num2str(jobnum) '-PSD-' num2str(ifonum)]);
close;

printmap(log10(DATA), psdjob.times, psdjob.f, 'Time','Frequency [Hz]', ['PSD ' num2str(ifonum)]);
colormap(parula), h2=colorbar;
title(['Map of PSD ' num2str(ifonum)]);
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_INIDVIDUAL-JOB-' num2str(jobnum) '-PSD-CUT-' num2str(ifonum)]);
close;

