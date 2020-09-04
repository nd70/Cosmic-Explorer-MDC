function params = pprocParamsCheck(params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that checks post processing parameters
% for radiometer analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that output directories are set
try 
   params.post_processing_outputFilePrefix;
   params.output_plot_dir_prefix;
catch
   error('Please specify where you want output to go.');
end

% check that absolute sigma cut is set
try 
   params.doAbsoluteSigmaCut;
catch
   params.doAbsoluteSigmaCut = 0;
end

% check that delta sigma cut is set
try 
   if ~params.doDeltaSigmaCut
   fprintf('WARNING: delta sigma cut turned off\n');
   end
catch
   fprintf('WARNING: delta sigma cut turned off\n');
   params.doDeltaSigmaCut = 0;
end

% check if badGPSTimesFile was passed
try 
   params.badGPSTimesFile;
catch
   params.badGPSTimesFile = '';
end
% is the data for narrowband radiometer or not?
try 
   params.doNarrowbandRadiometer;
catch
   params.doNarrowbandRadiometer = 0;
   fprintf('Note that Narrowband Radiometer is set to false!\n');
end


% check that skydirection name for file saving and 
% name for title of plots are set
if params.doNarrowbandRadiometer
   try
      params.skyDirectionName;
      params.skyDirectionTitle;
   catch
      error('set sky drection name and how you would like it displayed in titles');
   end
end
% save data for each epoch if data is broken into epochs?
try
   params.saveEpochData;
catch
   params.saveEpochData = 0;
end

% analysis file containing all params and jobs files to combine
try 
   params.analysisFile;
catch
   error('Analysis File not given.\n');
end

% interferometers used in analysis
try 
   params.ifos;
catch
   error('Interferometers in analysis are not set\n');
end

% amplitude scaling factors
try 
   params.amplitude_scalings;
   if ~(length(params.amplitude_scalings)==length(params.ifos))
      error('either too many or too few amplitude scaling factors given number of ifos entered\n');
   end
catch
   params.amplitude_scalings = ones(1,length(params.ifos));
   fprintf('Amplitude scaling factors not set...setting them to 1\n');
end

% calibration errors
try 
   params.errs;
   if ~(length(params.errs)==length(params.ifos))
      error('either too many or too few amplitude scaling factors given number of ifos entered\n');
   end
catch
   params.errs = zeros(1,length(params.ifos));
end

% confidence level
try
   params.conf;
catch
   params.conf = 0.90;
   fprintf('WARNING: confidence level not set. Setting it to .90\n');
end

% bias correction for using data from adjacent segments to estimate variance
try
   params.bias_factor;
   catch
   params.bias_factor = 1;
   fprintf('WARNING: bias correction for using data to estimate variance is not set. Setting it to 1\n');
end

%flows
if params.doNarrowbandRadiometer
   try
      params.flows;
      params.fhighs;
      params.deltaF;
      if ~(length(params.flows) == length(params.fhighs))
         error('length of flows and fhighs do not match\n');
      end
   catch
      error('flows or fhighs for deltaF not set\n');
   end
end
%doCondor
try
   params.doCondor;
catch
   params.doCondor = 0;
end

%jobsToCombine
try
   params.jobsToCombine;
catch
   fprintf('WARNING: params.jobsToCombine not set.\n Setting it to 250 stochastic jobs to combine\n');
   params.jobsToCombine = 250;
end

% where to find prior results?
try
   params.prior_results_file;
catch
  params.prior_results_file = '';
end

% make plots
try
   params.doPlots;
catch
   params.doPlots = 0;
end

try 
   params.dagFullName;
catch
   if params.doCondor
   fprintf('do Condor is on and you havent specified where dag should go. putting in this directory\n');
   params.dagFullName = './pproc_rad';
   end
end

try 
   params.exeDir;
catch
   if params.doCondor
      fprintf('do Condor is on and executable directory not set. setting it to this directory.');
      params.exeDir = './';
   end
end
try 
    params.segmentDuration;
catch
    error('segmentDuration variable not set in params file');
end

try
    params.skippedSkyDirections;
catch
    params.skippedSkyDirections = 0;
end

try
    params.doExtraNotching;
catch
    fprintf('doExtraNotching not set. Default is false\n');
    params.doExtraNotching = false;
end
