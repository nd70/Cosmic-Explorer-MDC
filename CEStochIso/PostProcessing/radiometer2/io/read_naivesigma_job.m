function [nsigma sigma] = read_naivesigma_file(params, jobnum);

if isstr(jobnum)
    jobnum = str2num(jobnum);
end
if isstr(params)
    params = readParamsFromFile(params);
end

[t, nsigma.data, sigma.data] = textread([params.outputFilePrefix '_naivesigmas.job' num2str(jobnum) '.trial1.dat'],'%f %f %f','commentstyle','matlab');
nsigma.times = t;
sigma.times = t;
