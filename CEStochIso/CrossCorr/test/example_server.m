stochastic_paths;
delete('example_server.out');
[success, msg, msgid] = rmdir('cachefiles_server', 's');
mkdir('cachefiles_server');
diary('example_server.out');
stochastic('example_server.txt', 'example_jobs.txt', 2);
diary off;
