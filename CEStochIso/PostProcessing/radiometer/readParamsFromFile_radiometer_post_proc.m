function params = readParamsFromFile_radiometer_post_proc(paramsFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  readParamsFromFile --- read in search parameters from a file
%
%  readParamsFromFile(paramsFile) reads in search parameters from a
%  file, returning the parameters in a structure.
%
%  Assumes parameters are given by name/value pair.
%
%  Routine written by Joseph D. Romano, John T. Whelan, Vuk Mandic.
%  Contact Joseph.Romano@astro.cf.ac.uk, john.whelan@ligo.org and/or
%  vmandic@ligo.caltech.edu
%
%  $Id: readParamsFromFile.m,v 1.38 2008-12-17 16:08:53 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% read in name/value pairs
[names,values] = ...
  textread(paramsFile, '%s %s\n', -1, 'commentstyle', 'matlab');

%% check that number of names and values are equal
if length(names)~=length(values)
  error('invalid parameter file');
end

%% loop over parameter names, assigning values to structure
for ii=1:length(names)

  switch names{ii}

    case 'doPlots'
      params.doPlots = str2num(values{ii});
  
    case 'post_processing_outputFilePrefix'
      params.post_processing_outputFilePrefix = values{ii};

    case 'output_plot_dir_prefix'
      params.output_plot_dir_prefix = values{ii};

    case 'skyDirectionName'
      params.skyDirectionName = regexp(values{ii},',','split');

    case 'skyDirectionTitle'
      params.skyDirectionTitle = regexp(values{ii},',','split');

    case 'skippedSkyDirections'
      params.skippedSkyDirections =transpose(str2num(values{ii})); 

    case 'doAbsoluteSigmaCut'
      params.doAbsoluteSigmaCut = str2num(values{ii});

    case 'doDeltaSigmaCut'
      params.doDeltaSigmaCut = str2num(values{ii});
      
    case 'deltaSigmaThreshold'
      params.deltaSigmaThreshold = str2num(values{ii});

    case 'numJobs'
      params.numJobs = str2num(values{ii});

    case 'ifos'
      params.ifos = regexp(values{ii},',','split');

    case 'epochs'
      params.epochs = values{ii};

    case 'deltaF'
      params.deltaF = str2num(values{ii});

    case 'doFullBandEpochPlots'
      params.doFullBandEpochPlots = str2num(values{ii});

    case 'doFullBandFullRunPlots'
      params.doFullBandFullRunPlots = str2num(values{ii});

    case 'flows'
      params.flows = transpose(str2num(values{ii}));

   case 'fhighs'
      params.fhighs = transpose(str2num(values{ii}));

    case 'scale_factor'
      params.scale_factor = str2num(values{ii});

    case 'analysisDir'
      params.analysisDir = values{ii};

    case 'exeDir'
      params.exeDir = values{ii};
      
    case 'dagNameFull'
      params.dagNameFull = values{ii};

    case 'analysisFile'
      params.analysisFile = values{ii};

    case 'badGPSTimesFile'
      params.badGPSTimesFile = values{ii};
      
    case 'doCondor'
      params.doCondor = str2num(values{ii});

   case 'dagFullName'
      params.dagFullName = values{ii};

    case 'jobsToCombine'
      params.jobsToCombine = str2num(values{ii});

    case 'errs'
      params.errs = transpose(str2num(values{ii}));

    case 'amplitude_scalings'
      params.amplitude_scalings = transpose(str2num(values{ii}));

    case 'conf'
      params.conf = str2num(values{ii});

    case 'plotTitlePrefix'
      params.plotTitlePrefix = values{ii};
      
    case 'bias_factor'
      params.bias_factor = str2num(values{ii});

    case 'saveEpochData'
      params.saveEpochData = str2num(values{ii});

    case 'Nra'
      params.Nra = str2num(values{ii});

    case 'Ndec1' 
      params.Ndec1 = str2num(values{ii});

    case 'doNarrowbandRadiometer'
      params.doNarrowbandRadiometer = str2num(values{ii});
   
   case 'numSkyDirections'
      params.numSkyDirections = str2num(values{ii});

   case 'matAvailable'
      params.matAvailable = str2num(values{ii});

   case 'pprocID'
      params.pprocID = values{ii};
  
   case 'segmentDuration'
      params.segmentDuration = str2num(values{ii});
   
   case 'doExtraNotching'
      params.doExtraNotching = str2num(values{ii});
   
   case 'pproc_notch_file'
      params.pproc_notch_file = values{ii};

   otherwise
      %% do nothing 

    end %% switch

end %% loop over parameter names

return
