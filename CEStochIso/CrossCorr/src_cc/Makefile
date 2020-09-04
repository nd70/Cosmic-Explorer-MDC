TRUNK = ../../../../trunk
UTIL = ../../../../../../admin/utilities

MCC = mcc \
  -I $(TRUNK)/Utilities \
  -I $(TRUNK)/Utilities/ligotools/matlab \
  -I $(UTIL)/Channel \
  -I $(UTIL)/detgeom/matlab \
  -I $(UTIL)/misc \
  -I $(UTIL)/not4toolbox/FTSeries \
  -R -nojvm -R -nodisplay -R -singleCompThread

# Including iofun causes compilation to give a segmentation fault but
# doesn't appear to be needed in recent versions of Matlab (like 2013a)
#  -I ${MATLAB_ROOT}/toolbox/matlab/iofun \

all : stochastic

# Since the build can succeed even with errors, it is safest
# to clean every time
stochastic : clean
	$(MCC) -m stochastic.m

clean :
	rm -f stochastic requiredMCRProducts.txt mccExcludedFiles.log run_stochastic.sh readme.txt
