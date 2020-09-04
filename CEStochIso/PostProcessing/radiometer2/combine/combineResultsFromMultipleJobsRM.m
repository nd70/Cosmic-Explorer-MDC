function multiplecombinedjobs = combineResultsFromMultipleJobs(pproc_params, jobs)

TAG = getFileTag(pproc_params, jobs);

if isstr(pproc_params)
    load(pproc_params);
end

try
    pproc_params.cut.flow;
catch
    pproc_params.cut.flow = 0;
end
try
    pproc_parmas.cut.fhigh;
catch
    pproc_parmas.cut.fhigh = 0;
end
SINGLEJOBS = {};
first = 1;
for ii = 1:length(jobs)
    % get combined data for a single job
    combinedjob = combineResultsFromSingleJobRM(pproc_params, jobs(ii));
    SINGLEJOBS{ii} = combinedjob;
    try
        if isnan(combinedjob)
            continue
        end
    catch
    end

    if first
        multiplecombinedjobs = combinedjob;
        first = 0;
        continue;
    end
    % combine data
    multiplecombinedjobs.pte.data = (multiplecombinedjobs.pte.data .* multiplecombinedjobs.sigma.data.^-2 + ...
                                    combinedjob.pte.data .* combinedjob.sigma.data.^-2);
    multiplecombinedjobs.sigma.data = (multiplecombinedjobs.sigma.data.^-2 + ...
                                       combinedjob.sigma.data.^-2) .^ (-0.5);
    multiplecombinedjobs.pte.data = multiplecombinedjobs.pte.data ./ multiplecombinedjobs.sigma.data.^-2;
    multiplecombinedjobs.sigma.times = [multiplecombinedjobs.sigma.times;combinedjob.sigma.times];
    multiplecombinedjobs.pte.times = [multiplecombinedjobs.pte.times;combinedjob.pte.times];
    multiplecombinedjobs.sigma.badtimes = [multiplecombinedjobs.sigma.badtimes;combinedjob.sigma.badtimes];
    multiplecombinedjobs.pte.badtimes = [multiplecombinedjobs.pte.badtimes;combinedjob.pte.badtimes];
    multiplecombinedjobs.pte.data(isnan(multiplecombinedjobs.pte.data)) = 0;
    multiplecombinedjobs.sigma.data(isnan(multiplecombinedjobs.sigma.data)) = inf;
end
% save data
if pproc_params.save_individual_jobs
    save([pproc_params.directory '/' pproc_params.prefix '_INDIVIDUAL-JOBS' TAG],'SINGLEJOBS','pproc_params','-v7.3');
end
if pproc_params.save_final_combined_jobs
    FINAL_COMBINED = multiplecombinedjobs
    save([pproc_params.directory '/' pproc_params.prefix '_COMBINED-JOBS' TAG],'FINAL_COMBINED','pproc_params','-v7.3');
end
