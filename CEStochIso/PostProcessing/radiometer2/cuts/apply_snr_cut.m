function [cutjob MASK] = apply_cut(job)

snr_norm = [];
mask1 = [];
for ii = 1:length(job.psd1.data(1,:))
    snr = job.pte.data(:,ii) ./ job.sigma.data(:,ii);
    norm1 = 1;%median(snr);
    snr_norm = [snr_norm (snr / norm1)];
    mask1 = [mask1 abs(snr) > 6];
end
MASK = mask1;
tot = sum(MASK(:));
perc = tot / length(snr_norm(:));
MASK = ~MASK;
fprintf('cut percentage: %f\n',perc * 100);
job.psd1.data = job.psd1.data .* MASK;
job.psd2.data = job.psd2.data .* MASK;
job.sigma.data = job.sigma.data .* MASK;
job.sigma.data(job.sigma.data == 0) = inf;
job.pte.data = job.pte.data .* MASK;
cutjob = job;
