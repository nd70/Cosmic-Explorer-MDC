function save_jobs_to_mat_file(pproc_params)

if isstr(pproc_params)
    load(pproc_params);
end


numFiles = ceil(pproc_params.numJobs / pproc_params.jobs_to_combine);
for i = 1:numFiles
    js = (i-1) * pproc_params.jobs_to_combine + [1:pproc_params.jobs_to_combine];
    if js(end) > pproc_params.numJobs
        js = js(find(js <= pproc_params.numJobs));
    end
    jobs = read_multiple_stoch_jobs(pproc_params, js);
    TAG = getFileTag(pproc_params, js);
    save([pproc_params.directory '/' pproc_params.prefix '_saved_stochastic_jobs' TAG], 'jobs', 'pproc_params', '-v7.3');
end
