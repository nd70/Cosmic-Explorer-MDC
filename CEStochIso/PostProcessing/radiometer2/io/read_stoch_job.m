%function [job] = read_stoch_job(prefix, flow, fhigh, deltaF, jobnum, skyDir)
function [job] = read_stoch_job(params, jobnum, skyDir)

if isstr(params)
    params = readParamsFromFile(params);
end


job.psd1 = read_psd_file(params, 1, jobnum);
job.psd2 = read_psd_file(params, 2, jobnum);
[job.pte job.sigma] = read_ccstats_job(params, jobnum, skyDir);
[job.nsigma job.theorsigma] = read_naivesigma_job(params, jobnum);
