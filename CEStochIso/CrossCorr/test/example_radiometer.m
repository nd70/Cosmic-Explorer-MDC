clear all;
stochastic_paths;

delete('example_radiometer.out');
diary('example_radiometer.out');
stochastic('example_radiometer_params.txt', [ 'example_jobs.txt' ], 3);
diary off;
