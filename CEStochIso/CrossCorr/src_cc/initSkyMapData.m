function initSkyMapData(injectTimeDomain,...
                        fsample1,fsample2,nResample1,nResample2,betaParam1,betaParam2,...
                        Hf,intLog,Lmax,ifo1,ifo2,gammaLM_coeffsPath,...
                        det1,det2,isSpH,coord1,coord2,map,flow,deltaF,numFreqs,...
                        t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
			t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
			bufforget,bufspline,MakeIncoherent)
%
%function initSkyMapData(injectTimeDomain,...
%                        fsample1,fsample2,nResample1,nResample2,betaParam1,betaParam2,...
%                        Hf,intLog,Lmax,ifo1,ifo2,gammaLM_coeffsPath,...
%                        det1,det2,isSpH,coord1,coord2,map,flow,deltaF,numFreqs,...
%                        t1,f1,R01,C01,alpha1,gamma1,ASQchannel1,alphaBetaFile1,calCavGainFile1,calResponseFile1,...
%                        t2,f2,R02,C02,alpha2,gamma2,ASQchannel2,alphaBetaFile2,calCavGainFile2,calResponseFile2,...
%                        bufforget,bufspline,MakeIncoherent)
%
%  initSkyMapData --- initializes the global variable SKY_MAP_MEMORY.
%                          this function needs to be called before getSkyMapData is called
%                          for the first time.
%
%  arguments: see below
%
%  SKY_MAP_MEMORY struct gets initialized with the following entries:
%             injectTimeDomain: if true, inject in time-domain; if false, inject in frequency domain
%             fsample1   - sample frequency of the time series 1
%             fsample2   - sample frequency of the time series 2
%             nResample1 - order of matlab resample routine for time series 1
%             nResample2 - order of matlab resample routine for time series 2
%             betaParam1 - beta parameter for matlab resample routine for time series 1
%             betaParam2 - beta parameter for matlab resample routine for time series 2
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             Lmax       - Lmax for simulated sky map if isSpH is true (may differ from the Lmax used for analysis)
%             ifo1, ifo2 - strings identifying the two detectors (e.g., H1, L1)
%             glm, g1lm, g2lm - spherical harmonic components of the overlap reduction functions (for detectors 12, 11, 22)
%             det1,det2  - detector structures containing position r and tensor d
%             isSpH      - Type of map: true: map contains complex spherical harmonics; false: map is pixel map
%             coord1     - either vector of l or right ascension in hours
%                          takes vector for multiple point sources, must have same length as coord2
%             coord2     - either vector of m or declination in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as coord1
%	      map        - map data; either complex spherical haromnics, or value at pixel
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
%             MakeIncoherent optional parameter to destroy coherence of point source
%                            0:  coherent point source
%                            1:  incoherent, but scaled with DC antenna acceptance
%                            2:  stationary noise, PowerSpec corresponds to the noise seen in both 
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%  $Id: initSkyMapData.m,v 1.8 2008-09-25 16:39:44 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  MakeIncoherent;
catch
  MakeIncoherent = 0;
end;

global SKY_MAP_MEMORY;

% initialize SKY_MAP_MEMORY
SKY_MAP_MEMORY.injectTimeDomain=injectTimeDomain;
SKY_MAP_MEMORY.fsample1=fsample1;
SKY_MAP_MEMORY.fsample2=fsample2;
SKY_MAP_MEMORY.nResample1=nResample1;
SKY_MAP_MEMORY.nResample2=nResample2;
SKY_MAP_MEMORY.betaParam1=betaParam1;
SKY_MAP_MEMORY.betaParam2=betaParam2;
SKY_MAP_MEMORY.Hf=Hf;
SKY_MAP_MEMORY.intLog=intLog;
SKY_MAP_MEMORY.Lmax=Lmax;
SKY_MAP_MEMORY.ifo1=ifo1;
SKY_MAP_MEMORY.ifo2=ifo2;
SKY_MAP_MEMORY.det1=det1;
SKY_MAP_MEMORY.det2=det2;
SKY_MAP_MEMORY.isSpH=isSpH;
SKY_MAP_MEMORY.coord1=coord1;
SKY_MAP_MEMORY.coord2=coord2;
SKY_MAP_MEMORY.map=map;
SKY_MAP_MEMORY.flow=flow;
SKY_MAP_MEMORY.deltaF=deltaF;
SKY_MAP_MEMORY.numFreqs=numFreqs;
SKY_MAP_MEMORY.t1=t1;
SKY_MAP_MEMORY.f1=f1;
SKY_MAP_MEMORY.R01=R01;
SKY_MAP_MEMORY.C01=C01;
SKY_MAP_MEMORY.alpha1=alpha1;
SKY_MAP_MEMORY.gamma1=gamma1;
SKY_MAP_MEMORY.ASQchannel1=ASQchannel1;
SKY_MAP_MEMORY.alphaBetaFile1=alphaBetaFile1;
SKY_MAP_MEMORY.calCavGainFile1=calCavGainFile1;
SKY_MAP_MEMORY.calResponseFile1=calResponseFile1;
SKY_MAP_MEMORY.t2=t2;
SKY_MAP_MEMORY.f2=f2;
SKY_MAP_MEMORY.R02=R02;
SKY_MAP_MEMORY.C02=C02;
SKY_MAP_MEMORY.alpha2=alpha2;
SKY_MAP_MEMORY.gamma2=gamma2;
SKY_MAP_MEMORY.ASQchannel2=ASQchannel2;
SKY_MAP_MEMORY.alphaBetaFile2=alphaBetaFile2;
SKY_MAP_MEMORY.calCavGainFile2=calCavGainFile2;
SKY_MAP_MEMORY.calResponseFile2=calResponseFile2;
SKY_MAP_MEMORY.bufstart=-1e20;
SKY_MAP_MEMORY.bufforget=bufforget;
SKY_MAP_MEMORY.bufdur=-bufspline;
SKY_MAP_MEMORY.bufspline=bufspline;
SKY_MAP_MEMORY.buf1=[];
SKY_MAP_MEMORY.buf2=[];
SKY_MAP_MEMORY.MakeIncoherent=MakeIncoherent;

% calculation of gammaLM's depend on whether it is a time-domain or freq-domain injection
% NOTE: we will be able to interpolate the gammaLM's later on if flow and numFreqs are 
% chosen so that flow=1/duration and fhigh is larger than that required by the time domain
% simulation
if injectTimeDomain == true

  duration = SKY_MAP_MEMORY.bufspline*2;
  deltaF = 1;
  N1 = fsample1/deltaF;
  N2 = fsample2/deltaF;
  if N1 >= N2 
    N = N1;
  else
    N = N2;
  end

  % discrete positive frequencies (not including DC component and Nyquist)
  if ( mod(N,2)== 0 )
    numFreqs = N/2+2;
  else
    numFreqs = (N-1)/2+3;
  end
  flow = 1/duration;

end

SKY_MAP_MEMORY.glm =calGammaLM(gammaLM_coeffsPath,[ifo1(1),ifo2(1)],Lmax,numFreqs,flow,deltaF);
SKY_MAP_MEMORY.g1lm=calGammaLM(gammaLM_coeffsPath,[ifo1(1),ifo1(1)],Lmax,numFreqs,flow,deltaF);
SKY_MAP_MEMORY.g2lm=calGammaLM(gammaLM_coeffsPath,[ifo2(1),ifo2(1)],Lmax,numFreqs,flow,deltaF);

return;
