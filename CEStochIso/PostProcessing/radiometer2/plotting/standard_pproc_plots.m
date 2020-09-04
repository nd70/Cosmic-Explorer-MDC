function standard_pproc_plots(pproc_params, job1, njobs)

if isstr(pproc_params)
    load(pproc_params);
end

LB=flipud(lbmap(256,'BrownBlue'));

TAG = getFileTag(pproc_params, [job1:job1+njobs-1])

load([pproc_params.directory '/' pproc_params.prefix '_COMBINED-JOBS' TAG]);



SNR = FINAL_COMBINED.pte.data./FINAL_COMBINED.sigma.data;
SNR(SNR==0) = NaN;
[n1,x1] = hist(SNR,20);

% plot snr spectrum
figure;
subplot(1,3,[1 2]);
plot(FINAL_COMBINED.pte.f,SNR);
xlabel('Frequency [Hz]');
ylabel('SNR');
title('SNR Spectrum')
subplot(1,3,3);
stairs(n1,x1);
xlabel('number');
ylabel('SNR');
title('SNR histogram');
print('-dpng',[pproc_params.directory '/' pproc_params.prefix '_SNR-SPECTROGRAM' TAG]);

fle = [pproc_params.directory '/' pproc_params.prefix '_loud_bins' TAG '.txt'];
fid = fopen(fle,'w+');
idxs = find(~isnan(SNR));
frequency = FINAL_COMBINED.pte.f(idxs);
SNR = SNR(idxs);
[SNR_sort, I] = sort(abs(SNR));
SNR2 = SNR(I);
frequency = FINAL_COMBINED.pte.f(I);

fprintf(fid, '# Loud Frequency Bins for Segment Duration: %d and Resolution: %4.4f and nbins: %d\n', pproc_params.segmentDuration, pproc_params.deltaF, numel(SNR_sort));
fprintf(fid, '# Frequency\tSNR\tSignificance\n');
for i = 0:50
    sig = 1-erf(SNR_sort(end-i)/sqrt(2)).^(numel(SNR_sort));
    fprintf(fid, '%4.4f\t\t%4.4f\t\t%4.4f\n',frequency(end-i), SNR2(end-i), sig);
end
fclose(fid)
