
paramsfile = '../../input/paramfiles/H1L1/H1L1_O1_150918_150926_297_1726_params.txt';

jobs_to_load = [1:163];
skyDirNum = 1;
jobs = read_multiple_stoch_jobs(paramsfile, jobs_to_load, skyDirNum);
save('test','jobs','-v7.3');
