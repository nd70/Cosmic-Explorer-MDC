function run_radiometer_postproc(post_proc_paramsFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper script for running stochastic radiometer
% post processing pipeline.
% It either loops over all jobs needed to run
% or it calls a perl script (makeDagForPostProc2.pl)
% that creates a .dag and .sub file, compiles the nece
% ssary code and then submits the dag and sub files
% to condor.
%
% written by patrick meyers. meyers@physics.umn.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in parameters
pproc_params = readParamsFromFile_radiometer_post_proc(post_proc_paramsFile);
%params       = readParamsFromFile(paramsFile);
pproc_params = pprocParamsCheck(pproc_params);

% read off top line of analysis summary file
[epoch_vec paramFiles jobFiles] = textread(pproc_params.analysisFile,'%s %s %s');
epochs = unique(epoch_vec);

% make sure output directory structure is set up. 
% extract directory structure from prefix
%create desired directory structure
[outputDir prefix ext] = fileparts(pproc_params.post_processing_outputFilePrefix);
[outputPlotDir prefix ext] = fileparts(pproc_params.output_plot_dir_prefix);
if ~exist(outputDir)
   [nouse1] = system(['mkdir -p ' outputDir]);
end
if ~exist(outputPlotDir);
   [nouse2] = system(['mkdir -p ' outputPlotDir]);
end
% check if the user wants to use condor or run 
% everything interactively in matlab
% WARNING: could take a long time if run interactively!
if pproc_params.doCondor
if ~pproc_params.matAvailable
   if ~exist(pproc_params.exeDir)
      [nouse] = system(['mkdir -p ' pproc_params.exeDir]);
   end

  % create command to run perl script that creates .dag and .sub files
  cmd2 = sprintf('./makeDagForPostProc2.pl %s %s %s %s %s',pproc_params.analysisFile,post_proc_paramsFile,num2str(pproc_params.jobsToCombine),pproc_params.dagFullName,pproc_params.exeDir);
  fprintf(cmd2);
  
  [junk2] = system(cmd2);


  [exeNAvailable,msg] = system(['ls ' pproc_params.exeDir 'radiometer_postproc']);
  if exeNAvailable
    fprintf('executable does not exist in this directory, attempting to compile code...\n');
    fprintf('trying to compile postprocessing code found here: \n');
    which radiometer_postproc;
    compile = sprintf('mcc -R -nodisplay -R -singleCompThread -m radiometer_postproc -d %s',pproc_params.exeDir);
    eval(compile);
    fprintf('done compiling code\n');
  end


  cmd3 = sprintf('condor_submit_dag -maxjobs 500 %s.dag',pproc_params.dagFullName);
  [junk3] = system(cmd3);
else % matAvailable
   if pproc_params.doNarrowbandRadiometer
      SumLocalRad(post_proc_paramsFile);
   else
      SumLocalRad_skymaps(post_proc_paramsFile);
   end % doNarrow 
end % matAvailable

else % do Condor
  % instead run interactively. 
  % we run over each individual parameter file.
if ~pproc_params.matAvailable
  for ll=1:length(paramFiles)
    % read jobs file. 
    [num startGPS endGPS dur] = textread(cell2mat(jobFiles(ll)),'%d %d %d %d');
    % use jobs file to figure out how many post processing jobs we need.
    for jj=1:ceil(length(startGPS)/pproc_params.jobsToCombine)
      % run radiometer_postproc for each of these jobs for this parameter file.
      paramStruct.paramsFile=cell2mat(paramFiles(ll));
      radiometer_postproc(paramStruct,post_proc_paramsFile,cell2mat(jobFiles(ll)),cell2mat(epoch_vec(ll)), jj);
    end % jj
  end % ll
end %mat available
  if pproc_params.doNarrowbandRadiometer
    SumLocalRad(post_proc_paramsFile);
  else
    SumLocalRad_skymaps(post_proc_paramsFile);
  end % doNarrow 
end% condor
