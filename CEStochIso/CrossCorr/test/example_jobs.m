stochastic_paths;
delete('example_jobs.out');
diary('example_jobs.out');
stochastic('example_params.txt', 'example_jobs.txt', 2);
diary off;
