%
% Run all jobs in the examnple job file to vaildate
% post-processing
%

stochastic_paths;
delete('example_server_all.out');
[success, msg, msgid] = rmdir('cachefiles_server', 's');
mkdir('cachefiles_server');
diary('example_server_all.out');
stochastic('example_server_all.txt', 'example_jobs.txt', 1);
stochastic('example_server_all.txt', 'example_jobs.txt', 2);
stochastic('example_server_all.txt', 'example_jobs.txt', 3);
diary off;
