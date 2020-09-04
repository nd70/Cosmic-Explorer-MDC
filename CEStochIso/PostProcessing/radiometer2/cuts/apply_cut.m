function [cutjob MASK baddsigtimes] = apply_cut(job, pproc_params)


try
    if pproc_params.cut.flow == 0
        pproc_params.cut.flow = job.pte.f(1);
    end
catch
    pproc_params.cut.flow = job.pte.f(1);
end

try
    if fhigh == 0
        pproc_params.cut.fhigh = job.pte.f(end);
    end
catch
    pproc_params.cut.fhigh = job.pte.f(end);
end
try
    pproc_params.doSkyMap;
catch
    pproc_params.doSkyMap = 0;
end

idx1 = find(job.pte.f == pproc_params.cut.flow);
idx2 = find(job.pte.f == pproc_params.cut.fhigh);

MASK = zeros(size(job.pte.data));
MASK_L = zeros(size(job.pte.data));
try 
    cutType = pproc_params.cut.type;
catch
    cutType = 'none'
end

if strcmp(cutType,'none')
    MASK = zeros(size(job.psd1.data));
elseif strcmp(cutType,'psd')
    psd1_norm = [];
    psd2_norm = [];
     for ii = 1:length(job.psd1.data(:,1))
        norm1 = median(job.psd1.data(ii,:));
        norm2 = median(job.psd2.data(ii,:));
        psd1_norm = [psd1_norm;job.psd1.data(ii,:) / norm1];
        psd2_norm = [psd2_norm;job.psd2.data(ii,:) / norm2];
        if and(ii >= idx1 , ii <= idx2)
            MASK(ii,:) = (job.psd1.data(ii,:) / norm1) > 7;
            MASK_L(ii,:) = (job.psd2.data(ii,:) / norm2) > 7;
            MASK = MASK .* MASK_L;
        end
    end
    A = MASK;
    tMASK = circshift(A,1);
    MASK = MASK+tMASK;
    tMASK = circshift(A,-1);
    MASK = MASK+tMASK;
    tMASK = circshift(A,2);
    MASK = MASK+tMASK;
    tMASK = circshift(A,-2);
    MASK = MASK+tMASK;
    tMASK = circshift(A,1,2);
    MASK = MASK+tMASK;
    tMASK = circshift(A,-1,2);
    MASK = MASK+tMASK;
    tMASK = circshift(A,2,2);
    MASK = MASK+tMASK;
    tMASK = circshift(A,-2,2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,1), 1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,1), -1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-1), -1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-1), 1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,2), 1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,2), -1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-2), -1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-2), 1, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,2), 2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,2), -2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-2), -2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-2), 2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,1), 2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,1), -2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-1), -2, 2);
    MASK = MASK+tMASK;
    tMASK = circshift(circshift(A,-1), 2, 2);
    MASK = MASK+tMASK;


    MASK(MASK>1) = 1;
    tot = sum(MASK(:));
    perc = tot / length(psd1_norm(:));
    %fprintf('cut percentage: %f\n',perc * 100);
elseif strcmp(cutType,'pte')
     pte_norm = [];
     for ii = 1:length(job.pte.data(1,:))
        dat = job.pte.data(~isnan(job.pte.data(:,ii)), ii);
        norm1 = median(abs(dat));
        pte_norm = [pte_norm abs(job.pte.data(:,ii)) / norm1];
        MASK = [MASK ((abs(job.pte.data(:,ii)) / norm1) > 7)];
    end
    tot = sum(MASK(:));
    perc = tot / length(pte_norm(:));
    %fprintf('cut percentage: %f\n',perc * 100);
elseif strcmp(cutType,'snr')
    snr_norm = [];
    for ii = 1:length(job.psd1.data(1,:))
        snr = job.pte.data(:,ii) ./ job.sigma.data(:,ii);
        snr_norm = [snr_norm snr];
        MASK = [MASK abs(snr) > 5]; % mask for bad snr pixels
    end
    tot = sum(MASK(:));
    perc = tot / length(snr_norm(:));
    %fprintf('cut percentage: %f\n',perc * 100);
else
    error('you specified a cut type ttat''s not available...\n');
end
if pproc_params.doSkyMap
    MASK = zeros(size(job.sigma.data));
end
MASK = ~MASK;

% delta sigma cut
MASK2 = and(job.theorsigma.data ./ job.nsigma.data < 1.2,...
            job.theorsigma.data ./ job.nsigma.data > 0.8);

%MASK2 = zeros(length(job.pte.times),1);
MASK2(1:length(job.theorsigma.times)) = abs(log10(job.theorsigma.data ./ job.nsigma.data)) < log10(1.2);
% [filetimes ifo fl fh] = textread('/home/meyers/sgwb/trunk/O1/radiometer/input/misc/badGPSTimes_a3_20-1726Hz.txt','%d %s %d %d');
% filetimes = filetimes(:,1);
% for t = 1:length(filetimes)
%     idx = find(job.theorsigma.times == filetimes(t));
%     if ~isempty(idx)
%         MASK2(idx) = 0;
%     end
% end

baddsigtimes = job.pte.times(find(~MASK2 .* job.pte.times));

% expand delta sigma mask to full band
[MASK2, Y] = meshgrid(MASK2, ones(length(MASK(:,1)),1)); 

% combine the two masks. 
% true is good segment, false is bad segment.
MASK = MASK .* MASK2;
if ~pproc_params.doSkyMap
    job.psd1.data = job.psd1.data .* MASK(:,1:length(job.psd1.data(1,:)));
    job.psd2.data = job.psd2.data .* MASK(:,1:length(job.psd1.data(1,:)));
end
job.sigma.data = job.sigma.data .* MASK(:,1:length(job.sigma.data(1,:)));
job.sigma.data(job.sigma.data == 0) = inf;
job.pte.data = job.pte.data .* MASK(:,1:length(job.pte.data(1,:)));
job.badtimes  = baddsigtimes;
cutjob = job;
