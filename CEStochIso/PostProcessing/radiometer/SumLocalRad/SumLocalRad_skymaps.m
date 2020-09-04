function SumLocalRad_skymaps(pproc_paramsFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pproc_params = readParamsFromFile_radiometer_post_proc(pproc_paramsFile);
pproc_params = pprocParamsCheck(pproc_params);
matPrefix = pproc_params.post_processing_outputFilePrefix;
temp = regexp(matPrefix,'/','split');

matDir = [];
for i = 1:length(temp)-1
matDir = [matDir temp{i} '/'];
end

%analysisFile = pproc_params.analysisFile;
%[epoch_vec,paramsFiles,jobsFiles] = textread(analysisFile,'%s %s %s');
mats_to_load = dir([matPrefix '*.mat']);

% [jobsFiles] = unique(jobsFiles);

pproc_params.Nra = 360;
pproc_params.Ndec1 = 181;
flows  = pproc_params.flows;
fhighs = pproc_params.fhighs;

E=[]; % initialize vector for 1/(sigma^2)
Y=[]; % initialize vector for ptEst/(sigma^2)
for ii=1:length(pproc_params.flows)
    E = [];
    Y = [];
  for jj=1:length(mats_to_load);
    pproc_dat = load([matDir mats_to_load(jj).name]);
%    try pproc_dat.flow; catch continue; end
%    try pproc_dat.paramsFile; catch continue; end
      
      % to add these two vectors together all post processing params must be the same between
      % file that's being loaded and those passed to this file
      % flow must match (fhigh whould be taken care of by the paramsFile check and
      % the pproc_params check)
  %    if (pproc_dat.flow == pproc_params.flows(ii) & any(strcmp(pproc_dat.paramsFile,paramsFiles)) & isequal(pproc_params,pproc_dat.pproc_params));
%	params = readParamsFromFile(pproc_dat.paramsFile);
%	scale_factor = determine_scale_factor(pproc_params.ifos,pproc_params.amplitude_scalings,[params.ifo1 params.ifo2]);
%
	try 
        pproc_dat.sky;
        catch
        continue;
        end
	if (isempty(E))
	  E = zeros(length(pproc_dat.sky(:,2)),1);
	  Y = zeros(length(pproc_dat.sky(:,2)),1);
	end % isempty E narrow
        try
	E = E + pproc_dat.sky(:,2).^-2;
	Y = Y + pproc_dat.sky(:,1).*pproc_dat.sky(:,2).^-2;
        catch
        end
  %    end % if flows match
    end % for jj
    sig = (E.^-0.5)*pproc_params.bias_factor;
    pte = Y./E;
    figure;
    histfit(pte./sig,50,'normal');
    xlabel('snr');
    ylabel('number of sky pixels');
    axis([-4 4 0 4000]);
    pretty;
    title(['SNR Histogram']);
    print('-depsc2',[pproc_params.output_plot_dir_prefix '_snr_histogram']);

    [h,p] = kstest(pte./sig);
    fprintf('%f is the p value for the ks test for the snr distribution\n',p);
    

    if ~isempty(E)
      map = [pte sig];
      plotMapAitoff_skymaps(map,pproc_params.Nra,pproc_params.Ndec1,pproc_params.output_plot_dir_prefix);
      pproc_params.final_data = true;
      save([pproc_params.post_processing_outputFilePrefix '_final_data'],'map','pproc_params','p');
    end
end % for ii
