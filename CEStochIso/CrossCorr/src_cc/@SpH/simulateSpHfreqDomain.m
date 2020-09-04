function [simCP, simP1, simP2] = simulateSpH(signalType, spectralAmplitude, ...
                                             spectralIndex, siderealTime, ...
                                             vSpH, doRandomize, modifyPSDs);

%
% simulates CSD and PSDs for point sources or multipole moments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set defaults for doRandomize, modifyPSD
try
  doRandomize;
catch
  doRandomize = false;
end;
try
  modifyPSDs;
catch
  modifyPSDs = false;
end;

% initialize plm's
lvec = vSpH.glm.lvec;
mvec = vSpH.glm.mvec;
plm = zeros(size(lvec));

% generate plms for desired signal type
switch signalType

  case 'pointSources'

    ra = [ 0 6 ];
    decl = [ 90 45 ];

    for ii=1:length(ra)
      pNP = deltaAtNorthPole(vSpH.Lmax);
      wRot = WignerRotation(vSpH.Lmax, (90-decl(ii))/180*pi);
      Z = zRotation(vSpH.Lmax, ra(ii)/12*pi);
      plm = plm + Z*wRot*pNP;
    end

  case 'monopole'
 
    l = 0;
    m = 0;

    ind1 = find(lvec==l);
    ind2 = find(mvec==m);
    ind = intersect(ind1,ind2);
    plm(ind) = 1;

  case 'dipole'
 
    %l = 1;
    %m = 0;
    %p = 1;

    l = 1;
    m = [1 -1];
    p = [-1/sqrt(2) 1/sqrt(2)];

    for ii=1:length(m)
      ind1 = find(lvec==l);
      ind2 = find(mvec==m(ii));
      ind = intersect(ind1,ind2);
      plm(ind) = p(ii);
    end

  case 'quadrupole'
 
    l = 2;
    m = 0;

    ind1 = find(lvec==l);
    ind2 = find(mvec==m);
    ind = intersect(ind1,ind2);
    plm(ind) = 1;
  
  otherwise
    error('unrecognized signal type');

end;

% randomize if desired
if doRandomize
  % modify plms by random plms
  random_plm = randPlm(vSpH.Lmax, true);
  plm = plm+100*random_plm;
end

% calculate simulated cross-power 
simCP = SpH2CrossPower(vSpH.glm, plm, siderealTime, ...
                         vSpH.numFreqs, vSpH.flow, vSpH.deltaF, ...
                         spectralAmplitude, spectralIndex);

% set simulated PSDs appropriately
if modifyPSDs
  % use simulated PSDs associated with plms
  simP1 = SpH2CrossPower(vSpH.g1lm, plm, siderealTime, ...
                           vSpH.numFreqs, vSpH.flow, vSpH.deltaF, ...
                           spectralAmplitude, spectralIndex);
  simP1 = abs(simP1); % make sure non-negative for randomized plm's
  simP2 = SpH2CrossPower(vSpH.g2lm, plm, siderealTime, ...
                           vSpH.numFreqs, vSpH.flow, vSpH.deltaF, ...
                           spectralAmplitude, spectralIndex);
  simP2 = abs(simP2); % make sure non-negative for randomized plm's

else
  % set to zeros
  simP1 = zeros(size(simCP));
  simP2 = zeros(size(simCP));
end

return

