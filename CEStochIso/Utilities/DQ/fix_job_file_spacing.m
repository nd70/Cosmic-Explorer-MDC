function fix_job_file(old_file)
% gwynne's code for fixing job file
%
%
data = dlmread(old_file);
outname = old_file;
fid=fopen(outname,'w');
for ii=1:length(data)
    spcs = '       ';
    nsdiff = data(ii,3)-data(ii,2);
    len = length(num2str(nsdiff));
    fprintf(fid, '   1  %i  %i%s%i\n', data(ii,2), data(ii,3), spcs(1:7-len), nsdiff);%% print a job file line in standard format accepted by stochastic.m
end
fclose(fid);
