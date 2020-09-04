function [pte, sigma] = read_ccstats_job(params, jobnum, skyDirection)
%function [pte, sigma] = read_ccstats_job(prefix, flow, fhigh, deltaF, jobnum, directionNum)

if isstr(params)
    params = readParamsFromFile(params);
end
if isstr(jobnum)
    jobnum = str2num(jobnum);
end

if isstr(skyDirection)
    skyDirection = str2num(skyDirection);
end

try 
    skyDirection;
catch
    skyDirection = 0;
end


T = []; 
sigma.data = []; 
pte.data = []; 
f = [params.flow:params.deltaF:params.fhigh];

for jj = 1:50
  file = [params.outputFilePrefix '_ccstatsSkySet' num2str(jj) '.job' num2str(jobnum) '.trial1.mat'];
  try 
    g = load(file);    
    g.Sky{1};
  catch 
    continue;
  end 
  if skyDirection==0
     i1 = 1;
     i2 = length(g.Sky{1}.data(:,2));
  else
     i1 = (skyDirection - 1) * length(f) + 1;
     i2 = i1 + length(f) - 1;
  end
  for kk = 1:length(g.Sky)
    T = [T;g.Sky{kk}.time];
    sigma.data = [sigma.data g.Sky{kk}.data(i1:i2,2)];
    pte.data = [pte.data g.Sky{kk}.data(i1:i2,1)];
  end 
end 
sigma.times = T;
pte.times = T;
sigma.f = f;
pte.f = f;
