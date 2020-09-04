
MPATH = $(shell which matlab | sed 's/bin.*$$//g' )

MCC = mcc -I ../Utilities -I ../../../../admin/utilities/Channel \
  -I ../lib/Clones/ \
  -I ../../../../admin/utilities/misc \
  -I ../../../../admin/utilities/not4toolbox/FTSeries \
  -I ../../../../admin/utilities/detgeom/matlab \
  -I ../../../../packages/burst/rStatistic \
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
src_cc/calCrossCorr.m \
src_cc/calOptimalFilter.m \
src_cc/constructFreqMask.m \
src_cc/normalizationIntegral.m \
src_cc/normalization.m \
src_cc/readCalibrationFromFiles.m \
src_cc/readParamsFromFile.m \
src_cc/readTimeSeriesData.m \
src_cc/readTSDataFromMatfile.m \
src_cc/stochastic.m \
src_cc/windowAndFFT.m \
src_cc/windowFactors.m \
../Utilities/calculateResponse.m \
../Utilities/coarseGrain.m \
../Utilities/convertResponse.m
	$(MCC) -m src_cc/stochastic.m

tidy :
	rm -f mccExcludedFiles.log *.c *.h

clean : tidy
	rm -f stochastic stochastic_main.c stochastic_mcc_component_data.c \
	stochastic.ctf; \
	rm -rf stochastic_mcr; \
	rm -rf bin

