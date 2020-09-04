% SET AUTOMATIC PARAMS FROM FILE % 

params = readParamsFromFile(pproc_params.paramsFile);

% parameters for params file % 
pproc_params.flow = params.flow;
pproc_params.segmentDuration = params.segmentDuration;
pproc_params.fhigh = params.fhigh;
pproc_params.deltaF = params.deltaF;
pproc_params.SkyPatternFile = params.SkyPatternFile;
pproc_params.outputFilePrefix = params.outputFilePrefix;
[a,b,c] = fileparts(pproc_params.outputFilePrefix);
pproc_params.prefix = b;
