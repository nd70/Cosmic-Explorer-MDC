MCC = mcc \
  -I ../CrossCorr \
  -I ../../../utilities/Clones \
  -I ../../../utilities/FTSeries \
  -I ${MATLAB_ROOT}/toolbox/matlab/graph2d \
  -I ${MATLAB_ROOT}/toolbox/matlab/iofun \
  -I ${MATLAB_ROOT}/toolbox/signal/signal \
  -I ${MATLAB_ROOT}/toolbox/matlab/specgraph \
  -R -nojvm \
  -m

exe : postProcessScriptFull combineResultsFromMultiplePairs
mex :
docs :
test_exe :

# Since the build can succeed even with errors, I'm just
# going to be lazy and clean every time
postProcessScriptFull : clean \
  postProcessScriptFull.m \
  runningPointEstimate.m \
  StatisticalAnalysisofResults_v2.m \
  FFTofPtEstIntegrand.m
	$(MCC) postProcessScriptFull.m

combineResultsFromMultiplePairs : clean combineResultsFromMultiplePairs.m
	$(MCC) combineResultsFromMultiplePairs.m

tidy :
	rm -f mccExcludedFiles.log *.c *.h

clean : tidy
	rm -f postProcessScriptFull postProcessScriptFull_main.c \
	postProcessScriptFull_mcc_component_data.c \
	rm -rf bin

