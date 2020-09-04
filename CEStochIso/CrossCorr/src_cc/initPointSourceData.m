function initPointSourceData(fsample,Hf,intLog,det1,det2,ra,decl,power,flow,deltaF,numFreqs,...
                             t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
			     t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
			     bufforget,bufspline,MakeIncoherent)
%  function initPointSourceData(fsample,Hf,intLog,det1,det2,ra,decl,power,flow,deltaF,numFreqs,...
%                               t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
%                               t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
%                               bufforget,bufspline,MakeIncoherent)
%
%  initPointSourceData --- initializes the global variable POINT_SOURCE_MEMORY.
%                          this function needs to be called before getPointSourceData is called
%                          for the first time.
%
%  arguments: see below
%
%  POINT_SOURCE_MEMORY struct gets initialized with the following entries:
%             fsample    - sample frequency of the time series
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             det1,det2  - detector structures containing position r and tensor d
%             ra         - right ascension in hours of source in the sky
%                          takes vector for multiple point sources, must have same length as decl
%             decl       - right ascension in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as ra
%	      power      - power emanating from corresponding point source (aka coefficient of Hf for a given source)
%			   takes vector for multiple point sources, must have same length as ra and decl
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
%                         from simulatePointSource
%             buf1        the actual buffer for IFO1
%             buf2        the actual buffer for IFO2
%             MakeIncoherent optional parameter to destroy coherence of point source
%                            0:  coherent point source
%                            1:  incoherent, but scaled with DC antenna acceptance
%                            2:  stationary noise, PowerSpec corresponds to the noise seen in both 
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  MakeIncoherent;
catch
  MakeIncoherent = 0;
end;

global POINT_SOURCE_MEMORY;

% initialize POINT_SOURCE_MEMORY
POINT_SOURCE_MEMORY.fsample=fsample;
POINT_SOURCE_MEMORY.Hf=Hf;
POINT_SOURCE_MEMORY.intLog=intLog;
POINT_SOURCE_MEMORY.det1=det1;
POINT_SOURCE_MEMORY.det2=det2;
POINT_SOURCE_MEMORY.ra=ra;
POINT_SOURCE_MEMORY.decl=decl;
POINT_SOURCE_MEMORY.power=power;
POINT_SOURCE_MEMORY.flow=flow;
POINT_SOURCE_MEMORY.deltaF=deltaF;
POINT_SOURCE_MEMORY.numFreqs=numFreqs;
POINT_SOURCE_MEMORY.t1=t1;
POINT_SOURCE_MEMORY.f1=f1;
POINT_SOURCE_MEMORY.R01=R01;
POINT_SOURCE_MEMORY.C01=C01;
POINT_SOURCE_MEMORY.alpha1=alpha1;
POINT_SOURCE_MEMORY.gamma1=gamma1;
POINT_SOURCE_MEMORY.ASQchannel1=ASQchannel1;
POINT_SOURCE_MEMORY.alphaBetaFile1=alphaBetaFile1;
POINT_SOURCE_MEMORY.calCavGainFile1=calCavGainFile1;
POINT_SOURCE_MEMORY.calResponseFile1=calResponseFile1;
POINT_SOURCE_MEMORY.t2=t2;
POINT_SOURCE_MEMORY.f2=f2;
POINT_SOURCE_MEMORY.R02=R02;
POINT_SOURCE_MEMORY.C02=C02;
POINT_SOURCE_MEMORY.alpha2=alpha2;
POINT_SOURCE_MEMORY.gamma2=gamma2;
POINT_SOURCE_MEMORY.ASQchannel2=ASQchannel2;
POINT_SOURCE_MEMORY.alphaBetaFile2=alphaBetaFile2;
POINT_SOURCE_MEMORY.calCavGainFile2=calCavGainFile2;
POINT_SOURCE_MEMORY.calResponseFile2=calResponseFile2;
POINT_SOURCE_MEMORY.bufstart=-1e20;
POINT_SOURCE_MEMORY.bufforget=bufforget;
POINT_SOURCE_MEMORY.bufdur=-bufspline;
POINT_SOURCE_MEMORY.bufspline=bufspline;
POINT_SOURCE_MEMORY.buf1=[];
POINT_SOURCE_MEMORY.buf2=[];
POINT_SOURCE_MEMORY.MakeIncoherent=MakeIncoherent;

return;
