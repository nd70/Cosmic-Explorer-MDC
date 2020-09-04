stochastic_paths;
delete('test_jobs.out');
diary('test_jobs.out');
stochastic('old/input/paramfiles/test_params.txt', 'old/input/jobfiles/test_jobs.txt', 1);
diary off;
