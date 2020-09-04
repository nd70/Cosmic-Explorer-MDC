function [ccSpectra, flow, deltaF, gpsTimes] = readCCspectraFromFile(filename)
%
%  readCCspectraFromFile -- read CC spectra from file
%
%  [ccSpectra, flow, deltaF, gpsTimes] = readCCspectraFromFile(filename)
%  returns a 2-d array of CC spectra read in from a file (spectra 
%  corresponding to different times are in different rows).  Also returned
%  are the GPS times and initial frequency and frequency spacing 
%  corresponding to the spectra.  The spectra are sorted according to GPS 
%  time.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: readCCspectraFromFile.m,v 1.1 2005-04-07 14:01:00 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % default return values
  ccSpectra = [];
  flow = [];
  deltaF = [];
  gpsTimes = [];
  
  checkFileExists('CC spectra file', filename);
  
  [pathstr, name, ext] = fileparts(filename);
  if (ext == '.mat');

    % Check if there are no intervals to include, which can happen if eg. the science segment
    % was too short. If there are no segments we return.
    load(filename, 'params');
    if (params.numIntervalsTotal < 1)
      return;
    end;

    load(filename, 'segmentStartTime', 'ccSpec');
    gpsTimes = segmentStartTime;
    % Frequency series characteristics are the same for all spectra so we use
    % the first one to set them.
    flow  = ccSpec(1, 1).flow;
    deltaF  = ccSpec(1, 1).deltaF;
    numFreqs = length(ccSpec(1, 1).data);

    numGPSTimes = length(gpsTimes);
    ccSpectra = zeros(numGPSTimes, numFreqs);
    for k = 1:numGPSTimes
      ccSpectra(k, :) = ccSpec(k, 1).data;
    end;
    clear segmentStartTime ccSpec;
  else
    % Read in data from file
    % For individual jobs files, there will usually be several spectra in a single file, grouped into
    % rows with the same GPS time.
    [gpsTimes, ccSigmas, freqs, spectra_real, spectra_imag] ...
        = textread(filename, '%f%f%f%f%f\n', -1, 'commentstyle', 'matlab');

    if isempty(gpsTimes)
      return;
    end
    
    % This pulls out the indices corresponding to each unique GPS time
    [gpsTimes, ind] = unique(gpsTimes);
    freqs = unique(freqs);
    
    numFreqs = length(freqs);
    flow = freqs(1);
    deltaF = freqs(2) - freqs(1);
    
    % create matrix with rows labeled by gps times and columns by frequencies
    ccSpectra = reshape(spectra_real+1i*spectra_imag, numFreqs, length(gpsTimes));
    ccSpectra = transpose(ccSpectra);
  end;

  % sort the arrays on gps start times (if not already sorted)
  [gpsTimes, ind] = sort(gpsTimes);
  ccSpectra = ccSpectra(ind,:);

return;

