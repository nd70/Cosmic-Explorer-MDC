function jobs = read_multiple_stoch_jobs(pproc_params, jobs_to_load)

if isstr(pproc_params)
    pproc_params = set_pproc_params(pproc_params);
end
try 
  pproc_params.skyDirection;
catch 
  pproc_params.skyDirection = 0;
end

params = readParamsFromFile(pproc_params.paramsFile);

for ii = 1:length(jobs_to_load)
   jobs{ii} = read_stoch_job(params, jobs_to_load(ii), pproc_params.skyDirection);
   jobs{ii}.jobnum = jobs_to_load(ii);
end
