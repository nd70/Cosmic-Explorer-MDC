function pproc_full(pproc_params)
if isstr(pproc_params)
    load(pproc_params);
end

% save_jobs_to_mat_file(pproc_params);
pproc_params.cut.type = 'none';
%    pproc_params.directory = '~/sgwb/trunk/O1/radiometer/output/NO-CUTS/';
pproc_params.skyDirection = 0;
mcj = combineResultsFromMultipleJobsRM(pproc_params,[1:pproc_params.numJobs]);
plot_skymaps(pproc_params);
%    plot_residuals_map(pproc_params, 1, pproc_params.numJobs);
%    standard_pproc_plots(pproc_params, 1, pproc_params.numJobs);


    % WITH CUTS % 
%     pproc_params.cut.type = 'psd';
%     pproc_params.directory = '~/sgwb/trunk/O1/radiometer/output/WITH-CUTS-WIDEST/';
%     %mcj = combineResultsFromMultipleJobsRM(pproc_params,[1:pproc_params.numJobs]);
%     plot_residuals_map(pproc_params, 1, pproc_params.numJobs);
%     standard_pproc_plots(pproc_params, 1, pproc_params.numJobs);

 exit;
