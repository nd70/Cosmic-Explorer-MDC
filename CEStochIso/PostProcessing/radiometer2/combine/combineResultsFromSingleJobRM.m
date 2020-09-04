function combinedjob = combineResultsFromSingleJobRM(pproc_params, jobnum)
%function combinedjob = combineResultsFromSingleJobRM(paramsfile, jobnum, directionNum, doCut, flow, fhigh)


try
    pproc_params.cut.fhigh;
catch
    pproc_params.cut.fhigh = 0;
end
try
    pproc_params.cut.flow;
catch
    pproc_params.cut.flow = 0;
end



%params = readParamsFromFile(paramsfile);

%JOB = read_stoch_job(params.outputFilePrefix,...
%                     params.flow, params.fhigh,...
%                     params.deltaF, jobnum, directionNum);
%JOB = read_stoch_job(pproc_params.outputfilePrefix,...
%                     pproc_params.flow, pproc_params.fhigh,...
%                     pproc_params.deltaF, jobnum,...
%                     pproc_params.skyDirection);
JOB = read_stoch_job(pproc_params.paramsFile, jobnum, pproc_params.skyDirection);
epsilon = 3 / 70;
% change to summable quantities
if isempty(JOB.sigma.data)
    combinedjob = NaN;
    return 
end
[JOB, MASK, badtimesdsig] = apply_cut(JOB, pproc_params);
% epsilon matrix
% circshift "BAD" matrix by 1
% to get segments where the previous segment was bad.
% Then add edge effects for first and last:


epsilon_t_minus_1 = 0.5 * epsilon * circshift(~MASK,-1,2);
epsilon_t_plus_1 = 0.5 * epsilon * circshift(~MASK,1,2);
epsilon_t_minus_1(:,1) = 0;
epsilon_t_plus_1(:,end) = 0;
epsilon_t_plus_1(:,1) = 0;
epsilon_t_minus_1(:,end) = 0; 
FISHER = JOB.sigma.data.^-2 - ...
         (JOB.sigma.data.^-2 + circshift(JOB.sigma.data.^-2,-1,2)).*epsilon_t_minus_1 - ...
         (JOB.sigma.data.^-2 + circshift(JOB.sigma.data.^-2,1,2)).*epsilon_t_plus_1;
X = JOB.pte.data.*JOB.sigma.data.^-2 - ...
         (JOB.sigma.data.^-2 + circshift(JOB.sigma.data.^-2,-1,2)).*epsilon_t_minus_1.*JOB.pte.data - ...
         (JOB.sigma.data.^-2 + circshift(JOB.sigma.data.^-2,1,2)).*epsilon_t_plus_1.*JOB.pte.data;

sigma = sum(FISHER,2).^-0.5;
pte = sum(X,2) ./ sigma.^-2;

pte(isnan(pte)) = 0;
sigma(isnan(sigma)) = inf;

combinedjob.pte.data = pte;
combinedjob.pte.f = JOB.pte.f;
combinedjob.pte.times = JOB.pte.times;
combinedjob.pte.badtimes = badtimesdsig;
combinedjob.sigma.data = sigma;
combinedjob.sigma.f = JOB.sigma.f;
combinedjob.sigma.times = JOB.pte.times;
combinedjob.sigma.badtimes = badtimesdsig;
