outfile = 'params/O1week1-15_20-297_100-40';
pproc_params.cut.type = 'none';
pproc_params.cut.flow = 1200;
pproc_params.cut.fhigh = 1300;
pproc_params.save_individual_jobs = true;
pproc_params.save_final_combined_jobs = true;
pproc_params.doSkyMap = false;
pproc_params.directory = '~/sgwb/trunk/O1/radiometer/output/weeks-1-15_100-40/';
pproc_params.skyDirection = 1;
pproc_params.jobs_to_combine = 100;
% pproc_params.paramsFile = '/home/meyers/sgwb/trunk/O1/radiometer/input/paramfiles/H1L1/H1L1_O1_150918_150926_20_297_params.txt';
pproc_params.paramsFile='/home/meyers/sgwb/trunk/O1/radiometer/input/paramfiles/H1L1/H1L1_O1_week1-15_20_297_100-40.txt';
pproc_params.jobsFile='/home/meyers/sgwb/trunk/O1/radiometer/input/jobfiles/JOB-FILE-1126623617-1135652417.txt';
%pproc_params.jobsFile = '~/sgwb/trunk/O1/radiometer/input/jobfiles/JOB-FILE-1126623617-1129383017.dat';
p = load(pproc_params.jobsFile);
pproc_params.numJobs = numel(p(:,1));
set_pproc_params
save(outfile, 'pproc_params');
