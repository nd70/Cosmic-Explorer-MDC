
MPATH = $(shell which matlab | sed 's/bin.*$$//g' )

BASEPATH  =../../../../../
STOCHPATH = ${BASEPATH}packages/stochastic/trunk/
CCPATH = ${BASEPATH}packages/stochastic/trunk/CrossCorr/

MCC = mcc -I ${STOCHPATH}Utilities -I ${BASEPATH}admin/utilities/Channel \
  -I ${STOCHPATH}lib/Clones/ \
  -I ${BASEPATH}admin/utilities/misc \
  -I ${BASEPATH}admin/utilities/not4toolbox/FTSeries \
  -I ${BASEPATH}admin/utilities/detgeom/matlab \
  -I ${BASEPATH}packages/burst/rStatistic \
  -I ${LIGOTOOLS}/matlab \
  -I $(MPATH)toolbox/signal/signal \
  -I $(MPATH)toolbox/matlab/iofun \
  -I $(MPATH)toolbox/matlab/graph2d \
  -R -nojvm

exe : stochastic
mex :
docs :
test_exe :

# Since the build can succeed even with errors, I'm just
# going to be lazy and clean every time
stochastic : clean \
${CCPATH}/src_cc/calCrossCorr.m \
${CCPATH}/src_cc/calOptimalFilter.m \
${CCPATH}/src_cc/constructFreqMask.m \
${CCPATH}/src_cc/normalizationIntegral.m \
${CCPATH}/src_cc/normalization.m \
${CCPATH}/src_cc/readCalibrationFromFiles.m \
${CCPATH}/src_cc/readParamsFromFile.m \
${CCPATH}/src_cc/readTimeSeriesData.m \
${CCPATH}/src_cc/readTSDataFromMatfile.m \
${CCPATH}/src_cc/stochastic.m \
${CCPATH}/src_cc/windowAndFFT.m \
${CCPATH}/src_cc/windowFactors.m \
${STOCHPATH}Utilities/calculateResponse.m \
${STOCHPATH}Utilities/coarseGrain.m \
${STOCHPATH}Utilities/convertResponse.m
	$(MCC) -m ${CCPATH}/src_cc/stochastic.m

tidy :
	rm -f mccExcludedFiles.log *.c *.h

clean : tidy
	rm -f stochastic stochastic_main.c stochastic_mcc_component_data.c \
	stochastic.ctf; \
	rm -rf stochastic_mcr; \
	rm -rf bin

