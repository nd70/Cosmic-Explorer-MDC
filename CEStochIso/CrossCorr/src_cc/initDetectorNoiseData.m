function initDetectorNoiseData(fsample1,fsample2,nResample1,nResample2,betaParam1,betaParam2,...
                        P1,P2,intLog,...
                        flow,deltaF,numFreqs,...
                        t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
			t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
			bufforget,bufspline)
%
%function initDetectorNoiseData(fsample1,fsample2,nResample1,nResample2,betaParam1,betaParam2,...
%                        P1,P2,intLog,...
%                        flow,deltaF,numFreqs,...
%                        t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
%    			 t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
%			 bufforget,bufspline)
%
%  initDetectorNoiseData --- initializes the global variable DETECTOR_NOISE_MEMORY.
%                            this function needs to be called before getDetectorNoiseData is called
%                            for the first time.
%
%  arguments: see below
%
%  DETECTOR_NOISE_MEMORY struct gets initialized with the following entries:
%             fsample1   - sample frequency of the time series 1
%             fsample2   - sample frequency of the time series 2
%             nResample1 - order of matlab resample routine for time series 1
%             nResample2 - order of matlab resample routine for time series 2
%             betaParam1 - beta parameter for matlab resample routine for time series 1
%             betaParam2 - beta parameter for matlab resample routine for time series 2
%             P1,2       - power spectra (one-sided, strain^2/Hz) in detectors 1,2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate P1,2 logarithmicly
%             flow             \
%             deltaF            > frequency information used for the main analysis
%             numFreqs         |
%             t1,f1,           \
%             R01, C01,         \
%             alpha1,            \
%             gamma1,             \
%             ASQchannel1          \
%             alphaBetaFile1        \
%             calCavGainFile1        \
%             calResponseFile1        \ Calibration information
%             t2,f2,                  |
%             R02, C02,              |
%             alpha2,               |
%             gamma2,              |
%             ASQchannel2         |
%             alphaBetaFile2     |
%             calCavGainFile2   |
%             calResponseFile2 |
%             bufstart    GPS start time of buffer
%             bufforget   memory duration of the buffer - this determines how far back
%                         the data can be recalled, i.e. the identical random series is produced
%             bufdur      number of seconds of data in the buffer - ready to be read out, <=bufforget
%             bufspline   duration in sec of 1/4 period of a cos/sin used to spline data
%                         2 x bufspline is also the time interal used to get new data
%                         from simulateSkyMap
%             buf1        the actual buffer for IFO1
%             buf2        the actual buffer for IFO2
%
%  Routine written by Stefan Ballmer, Joe Romano
%  Contact sballmer@ligo.mit.edu
%
%  $Id: initDetectorNoiseData.m,v 1.2 2008-09-25 16:39:44 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DETECTOR_NOISE_MEMORY;

% initialize DETECTOR_NOISE_MEMORY
DETECTOR_NOISE_MEMORY.fsample1=fsample1;
DETECTOR_NOISE_MEMORY.fsample2=fsample2;
DETECTOR_NOISE_MEMORY.nResample1=nResample1;
DETECTOR_NOISE_MEMORY.nResample2=nResample2;
DETECTOR_NOISE_MEMORY.betaParam1=betaParam1;
DETECTOR_NOISE_MEMORY.betaParam2=betaParam2;
DETECTOR_NOISE_MEMORY.P1=P1;
DETECTOR_NOISE_MEMORY.P2=P2;
DETECTOR_NOISE_MEMORY.intLog=intLog;
DETECTOR_NOISE_MEMORY.flow=flow;
DETECTOR_NOISE_MEMORY.deltaF=deltaF;
DETECTOR_NOISE_MEMORY.numFreqs=numFreqs;
DETECTOR_NOISE_MEMORY.t1=t1;
DETECTOR_NOISE_MEMORY.f1=f1;
DETECTOR_NOISE_MEMORY.R01=R01;
DETECTOR_NOISE_MEMORY.C01=C01;
DETECTOR_NOISE_MEMORY.alpha1=alpha1;
DETECTOR_NOISE_MEMORY.gamma1=gamma1;
DETECTOR_NOISE_MEMORY.ASQchannel1=ASQchannel1;
DETECTOR_NOISE_MEMORY.alphaBetaFile1=alphaBetaFile1;
DETECTOR_NOISE_MEMORY.calCavGainFile1=calCavGainFile1;
DETECTOR_NOISE_MEMORY.calResponseFile1=calResponseFile1;
DETECTOR_NOISE_MEMORY.t2=t2;
DETECTOR_NOISE_MEMORY.f2=f2;
DETECTOR_NOISE_MEMORY.R02=R02;
DETECTOR_NOISE_MEMORY.C02=C02;
DETECTOR_NOISE_MEMORY.alpha2=alpha2;
DETECTOR_NOISE_MEMORY.gamma2=gamma2;
DETECTOR_NOISE_MEMORY.ASQchannel2=ASQchannel2;
DETECTOR_NOISE_MEMORY.alphaBetaFile2=alphaBetaFile2;
DETECTOR_NOISE_MEMORY.calCavGainFile2=calCavGainFile2;
DETECTOR_NOISE_MEMORY.calResponseFile2=calResponseFile2;
DETECTOR_NOISE_MEMORY.bufstart=-1e20;
DETECTOR_NOISE_MEMORY.bufforget=bufforget;
DETECTOR_NOISE_MEMORY.bufdur=-bufspline;
DETECTOR_NOISE_MEMORY.bufspline=bufspline;
DETECTOR_NOISE_MEMORY.buf1=[];
DETECTOR_NOISE_MEMORY.buf2=[];

return;
