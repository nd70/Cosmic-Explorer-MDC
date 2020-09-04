function SumLocalRad(pproc_paramsFile)
% -------------------------------------------------------
% SumLocalRad(pproc_paramsFile)
% Function that collects output from condor jobs
% combines them into one final .mat file
% creates diagnostic plots (by calling plot_rad.m)
% -------------------------------------------------------
%
% Parameters:
% -----------
%   pproc_paramsFile : str
%       file that contains post processing parameters
%
% Returns:
% -------
%   NONE : nothing output directly. output saved to files
%       specified pproc_paramsFile
%
% CONTACT : patrick.meyers@ligo.org
% -------------------------------------------------------

% read and check parameter files
pproc_params = readParamsFromFile_radiometer_post_proc(pproc_paramsFile);
pproc_params = pprocParamsCheck(pproc_params);

% break down directory structure
% search out mat files
matPrefix = pproc_params.post_processing_outputFilePrefix;
temp = regexp(matPrefix,'/','split');
matDir = [];
for i = 1:length(temp)-1
    matDir = [matDir temp{i} '/'];
end
mats_to_load = dir([matPrefix '*.mat']);

% read in analysis file
analysisFile = pproc_params.analysisFile;
[epoch_vec,paramsFiles,jobsFiles] = textread(analysisFile,'%s %s %s');
% [jobsFiles] = unique(jobsFiles);

flows  = pproc_params.flows;
fhighs = pproc_params.fhighs;
deltaF = pproc_params.deltaF;
ff = [];

% set up full vector of frequencies
% used in total analysis
for ii = 1:length(flows)
    if ii > 1
        if flows(ii) == fhighs(ii-1)
            ff = [ff (flows(ii) + deltaF):deltaF:fhighs(ii)];
        else
            ff = [ff flows(ii):deltaF:fhighs(ii)];
        end % flows ==
    else
        ff = [ff flows(ii):deltaF:fhighs(ii)];
    end% ii>1
end% for ...
nbins = length(ff);

% number of frequencies in EACH INDIVIDUAL BAND
% used in analysis
nfreqs = []; % number we WANT from each range
nfreqs_accounting = [];% number EXPECT from each range
for ii = 1:length(flows)
    if ii < length(flows) & flows(ii+1) == fhighs(ii)
        % if there's overlap between frequency arrays from
        % two different runs, cut off end of lower range
        nfreqs = [nfreqs length(flows(ii):deltaF:(fhighs(ii)-deltaF))];
        nfreqs_accounting = [nfreqs_accounting length(flows(ii):deltaF:fhighs(ii))];
    else
        nfreqs = [nfreqs length(flows(ii):deltaF:fhighs(ii))];
        nfreqs_accounting = [nfreqs_accounting length(flows(ii):deltaF:fhighs(ii))];
    end % if ii ...
end% for ii ...

E=[]; % initialize vector for 1/(sigma^2)
Y=[]; % initialize vector for ptEst/(sigma^2)
for mm=1:pproc_params.numSkyDirections
    if any(pproc_params.skippedSkyDirections == mm)
        continue;
    end
        E_dir=[]; % initialize vector for 1/(sigma^2)
        Y_dir=[]; % initialize vector for ptEst/(sigma^2)
    for ii=1:length(pproc_params.flows)
        E_narrow = [];
        Y_narrow = [];
        for jj=1:length(mats_to_load);
            pproc_dat = load([matDir mats_to_load(jj).name]);
            try pproc_dat.flow; catch continue; end

            try pproc_dat.paramsFile; catch continue; end
      
      % to add these two vectors together all post processing params must be the same between
      % file that's being loaded and those passed to this file
      % flow must match (fhigh whould be taken care of by the paramsFile check and
      % the pproc_params check)
            if (pproc_dat.flow == pproc_params.flows(ii) & any(strcmp(pproc_dat.paramsFile,paramsFiles)) & isequal(pproc_params.pprocID,pproc_dat.pproc_params.pprocID));

                pproc_dat.pproc_params.deltaF = pproc_params.deltaF;

                % read in params from the parameter file used in this analysis
        	params = readParamsFromFile(pproc_dat.paramsFile);
                
                % scale factor needed based on calibration doc
        	scale_factor = determine_scale_factor(pproc_params.ifos,pproc_params.amplitude_scalings,[params.ifo1 params.ifo2]);

        	if (isempty(E_narrow))
                    E_narrow = zeros(length(pproc_dat.sky(:,2)),1);
                    Y_narrow = zeros(length(pproc_dat.sky(:,2)),1);
        	end % isempty E narrow
                % start index of what we want to pull out of data.
                % this requires using frequency array that was saved
                startindx = 1+(mm-1)*nfreqs_accounting(ii);
                if pproc_params.doExtraNotching
               	    [E_add Y_add] = dataReadOut(pproc_dat,flows(ii),fhighs(ii),pproc_dat.pproc_params,pproc_params.pproc_notch_file,startindx);
                else
               	    [E_add Y_add] = dataReadOut(pproc_dat,flows(ii),fhighs(ii),pproc_dat.pproc_params,'None',startindx);
                end
        	E_add(isnan(E_add)) = 0;
        	Y_add(isnan(Y_add)) = 0;

                % end index is nfreqs away (the frequencies we WANT to pull out)
                endindx = (startindx-1)+nfreqs(ii);
        	E_narrow = E_narrow(1:nfreqs(ii)) + E_add(1:nfreqs(ii));
        	Y_narrow = Y_narrow(1:nfreqs(ii)) + Y_add(1:nfreqs(ii));

            end % if flows match
        end % for jj
        if ~isempty(E_narrow)
            E_dir = [E_dir;E_narrow(1:nfreqs(ii))];
            %E_dir = [E_dir;E_narrow(startindx:endindx-remove)];
            Y_dir = [Y_dir;Y_narrow(1:nfreqs(ii))];
            %Y_dir = [Y_dir;Y_narrow(startindx:endindx-remove)];
        end % ~isempty
  end % for ii
  E = [E E_dir];
  Y = [Y Y_dir];
end
    if ~(length(E(:,1)) == nbins) || ~(length(Y(:,1)) == nbins);
        error('ERROR: length of data doesnt match number of frequency bins');

    end % length check
    plot_rad(pproc_params,Y,E,ff,scale_factor);
end
%end % loop over mm

function [E Y] = dataReadOut(data,flow,fhigh,pproc_params,notchfile,startindx);
% dat is two columns. first is Y second
% is sigma
%
% H1L1_masterA.txt
% then 12 is index of the A.
% need epoch for pulsar notches...
ff2=flow:data.pproc_params.deltaF:fhigh;
ff = [];
for ii = 1:(data.pproc_params.numSkyDirections - length(pproc_params.skippedSkyDirections))
    ff = [ff ff2];
end

% if no notch file, then freq mask will
% be set to zero
if strcmp(notchfile,'None')
    extra = [0 0];
else
    extra = load(notchfile);
end

% first column center bin
% other second column bins to remove
f_center = [extra(:,1)];
nBinsToRemove = [extra(:,2)];

fmask = [];
% set up new frequency mask based on
% cuts added
for ii=1:length(f_center);
    if (mod(nBinsToRemove(ii),2) ==0)
        error('Extra notching doesn''t handle even number of bins\n');
    end
    cut_width = (nBinsToRemove(ii) - 1)/2;
    cut_low = f_center(ii) - cut_width*data.pproc_params.deltaF - data.pproc_params.deltaF / 2;
    cut_high = f_center(ii) + cut_width*data.pproc_params.deltaF + data.pproc_params.deltaF / 2;
    fmask = [fmask cut_low:data.pproc_params.deltaF:cut_high];
end % for

% read out data and apply mask
if ~isempty(data.sky)
    for ii = 1:length(ff2)
        idx = ii-1+startindx;
        if any(fmask == ff2(ii))
            Y(ii) = NaN;
            sig(ii) = NaN;
        else
            Y(ii) = data.sky(idx,1);
            sig(ii) = data.sky(idx,2);
        end % if
    end % for
    Y = (Y./(sig.^2))';
    E = (sig.^(-2))';
else 
    fprintf('%d\n',length(ff2));
    Y=nan(length(ff2)*data.pproc_params.numSkyDirections+data.pproc_params.numSkyDirections,1);
    E=nan(length(ff2)*data.pproc_params.numSkyDirections+data.pproc_params.numSkyDirections,1);
end
end %fxn
