function map=ccStatReadout(det1, det2, GPS, source, ccStat, ccVar, tauUnphys)
%  function map=ccStatReadout(det1, det2, GPS, source, ccStat, ccVar)
%
%  reads out the propper ccStat and ccVar value for each source
%  and scales it with gamma0
%
%  Arguments: det1,det2  -  structures for the 2 detectors
%             GPS        -  GPS time used to calculate antenna pattern
%             source     -  Nx2 matrix containing right ascension and declination
%                           for all N points to be looked at.
%             ccStat     -  time series of ccStats vs time shift
%             ccVar      -  corresponding sigma^2
%             tauUnphys  -  Additional constant time shift to hide real result
%
%  Output:    map        -  Nx2 matrix containing ccStat and sigma for each source
%
%  Note:      This routine was optimized for speed. The naive approach of
%             calling orfIntegrandSymbolic N times was 100 (!) times
%             slower than this. The average execution time per sky position
%             is now down to 12 microseconds (on my laptop...).
%
%
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try 
  tauUnphys;
catch
  tauUnphys = 0;
end;

M=size(source,1);
map=zeros(M,2);
time=GPStoGreenwichMeanSiderealTime(GPS);
sigma=sqrt(ccVar);
N=length(ccStat.data);
s = det2.r - det1.r;
% distance between sites
distance = norm(s);
% unit vector Omega
c=299792458;
w=pi/12;
for ii=1:M
  psi=w*(time-source(ii,1));
  theta=-pi/2+pi/180*source(ii,2);
  ctheta=cos(theta);                stheta=sin(theta);
  cpsi  =cos(psi);                  spsi  =sin(psi);
  Omega1=-cpsi*stheta;              Omega2=spsi*stheta;                 Omega3=ctheta;
  pp11=ctheta^2*cpsi^2-spsi^2;      pp12=-(ctheta^2+1)*cpsi*spsi;       pp13=ctheta*cpsi*stheta;
                                    pp22=ctheta^2*spsi^2-cpsi^2;        pp23=-ctheta*spsi*stheta;
                                                                        pp33=stheta^2;
  px11=-2*ctheta*cpsi*spsi;         px12=ctheta*(spsi^2-cpsi^2);        px13=-stheta*spsi;
                                    px22=2*ctheta*cpsi*spsi;            px23=-stheta*cpsi;
                                                                       %px33=0;
  g.F1p=det1.d(1,1)*pp11+2*det1.d(1,2)*pp12+2*det1.d(1,3)*pp13+det1.d(2,2)*pp22+2*det1.d(2,3)*pp23+det1.d(3,3)*pp33;
  g.F2p=det2.d(1,1)*pp11+2*det2.d(1,2)*pp12+2*det2.d(1,3)*pp13+det2.d(2,2)*pp22+2*det2.d(2,3)*pp23+det2.d(3,3)*pp33;
  g.F1x=det1.d(1,1)*px11+2*det1.d(1,2)*px12+2*det1.d(1,3)*px13+det1.d(2,2)*px22+2*det1.d(2,3)*px23;%+det1.d(3,3)*px33;
  g.F2x=det2.d(1,1)*px11+2*det2.d(1,2)*px12+2*det2.d(1,3)*px13+det2.d(2,2)*px22+2*det2.d(2,3)*px23;%+det2.d(3,3)*px33;
  % divide by 2 since Hf is defined as sum over pol of power
  g.gamma0=(g.F1p * g.F2p + g.F1x * g.F2x)/2;	     
  g.tau=(Omega1*s(1)+Omega2*s(2)+Omega3*s(3))/c;
  % need a minus sign for the readout
  n=1+N/2 - (g.tau+tauUnphys)/ccStat.deltaT;
  % do linear interpolation
  np=ceil(n); nm=floor(n);
  wp=n-nm;    wm=1-wp;
  map(ii,1)=( wm*ccStat.data(nm) + wp*ccStat.data(np) )/g.gamma0;
  map(ii,2)=sigma/abs(g.gamma0);

  % lines for testing against slower code
  %gt=orfIntegrandSymbolic(det1,det2,time,source(ii,1),source(ii,2));
  %v=[norm(g.tau-gt.tau),norm(g.F1p-gt.F1p),norm(g.F1x-gt.F1x),norm(g.F2p-gt.F2p),norm(g.F2x-gt.F2x)];
  %fprintf('%d\t%d\t%d\t%d\t%d\n',v);
end

