function vSpH=dumpSpH(vSpH)

% dump all fields of SpH to screen



  fprintf('detectorPair:       %s\n',vSpH.detectorPair);
  fprintf('Lmax:               %d\n',vSpH.Lmax);
  fprintf('numFreqs:           %d\n',vSpH.numFreqs);
  fprintf('flow:               %f\n',vSpH.flow);
  fprintf('deltaF:             %f\n',vSpH.deltaF);
%  fprintf('H:                  %s\n',vSpH.H);
%  fprintf('mask:               %s\n',vSpH.mask);
  fprintf('w1w2bar:            %f\n',vSpH.w1w2bar);
  fprintf('w1w2squaredbar:     %f\n',vSpH.w1w2squaredbar);
  fprintf('segmentDuration:    %f\n',vSpH.segmentDuration);
  fprintf('FreqIntFlag:        %d\n',vSpH.FreqIntFlag);
  fprintf('prefix:             %s\n',vSpH.prefix);
  fprintf('suffix:             %s\n',vSpH.suffix);
  fprintf('maxSegmentsPerFile: %d\n',vSpH.maxSegmentsPerFile);
  fprintf('setPrefix:          %s\n',vSpH.setPrefix);
  fprintf('setNumber:          %d\n',vSpH.setNumber);
  fprintf('setOffset:          %f\n',vSpH.setOffset);
  fprintf('currentFilename:    %s\n',vSpH.currentFilename);
%  fprintf('out:                %s\n',vSpH.out);
%  fprintf('glm:                %s\n',vSpH.glm);
vSpH.out
