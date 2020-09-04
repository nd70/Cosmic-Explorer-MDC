function val = getFileTag(pproc_params, jobs)

if isstr(pproc_params)
    load(pproc_params);
end
try
    jobs;
catch
    jobs = [1:pproc_params.numJobs];
end

if pproc_params.doSkyMap
    pproc_params.skyDirection = 0;
end
val = ['-' num2str(pproc_params.skyDirection) '-' num2str(jobs(1)) '-' num2str(jobs(end))];
