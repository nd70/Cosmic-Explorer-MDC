function [psdjob] = read_psd_file(params, ifonum, jobnum)

if isstr(ifonum)
    ifonum = str2num(ifonum);
end
if isstr(jobnum)
    jobnum = str2num(jobnum);
end
if isstr(params)    
    params = readParamsFromFile(params);
end

file = [params.outputFilePrefix '_psd' num2str(ifonum) '.job' num2str(jobnum) '.trial1.dat'];
try
[times sigs freqs psds] = textread(file,'%f %f %f %f','commentstyle','matlab');
catch
keyboard
end

flow_indexes = find(freqs == params.flow);
fhigh_indexes = find(freqs == params.fhigh);

psdjob.data = [];
psdjob.times = [];
psdjob.f = [params.flow:params.deltaF:params.fhigh];

for i = 1:length(flow_indexes)
    i1 = flow_indexes(i);
    i2 = fhigh_indexes(i);
    psdjob.data = [psdjob.data psds(i1:i2)];
    if times(i1) ~= times(i2)
        error('something is wrong, times don''t match within frequency band');
    end
    psdjob.times = [psdjob.times times(i1)];
    psdjob.f = freqs(i1:i2);
end
