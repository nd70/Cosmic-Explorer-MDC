function [gpsTimes, naiveSigmas, ccSigmas] = readSigmasFromFile(fileList)
%
%  readSigmasFromFile -- read naive sigmas from file
%
%  [ gpsTimes, naiveSigmas, ccSigmas ] = readSigmasFromFile(fileList)
%
%  returns column vectors containing GPS times, naive sigmas and CC sigmas
%  from a file or cell array of files. The values are sorted according to GPS time.
%
%  The argument fileList may represent a single file name or a list of filenames
%  contained in a cell array eg. fileList could be
%
%     'filename.dat'
%     { 'filename.dat' }
%     { 'filename1.dat', 'filename2.dat', ... }
%
% The files may be text files with the format
%
%    <gpsTime> <naiveSigma> <ccSigma>
%
%  or Matlab files produced by stochastic.m using the writeOutputToMatFile flag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gpsTimes = [];
  naiveSigmas = [];
  ccSigmas = [];

  % If the file list is just a single filename, convert
  % it to a cellarray with one element.
  if (~iscell(fileList))
    fileList = { fileList };
  end;

  for f = fileList
    f = f{1};
    checkFileExists('Naive sigma file', f);

    [pathstr, name, ext] = fileparts(f);
    if (ext == '.mat');

      % Check if there are no intervals to include, which can happen if eg. the science segment
      % was too short. If there are no segments we continue to the next file
      load(f, 'params');
      if (params.numIntervalsTotal < 1)
        continue;
      end;

      load(f, 'segmentStartTime', 'naiSigma', 'ccSigma');
      gTmp = segmentStartTime;
      % Stats are stored in the .mat file as ccSigma(<segment index>, <trial index>). Post processing
      % doesn't handle multiple trials so we always use trial 1.
      nsTmp = naiSigma(:, 1);
      ccsTmp = ccSigma(:, 1);
      clear segmentStartTime naiSigma ccSigma;
    else
      % Read from text file
      [gTmp, nsTmp, ccsTmp] = textread(f, '%f%f%f\n', -1, 'commentstyle', 'matlab');

      % Check for empty results, which can happen if eg. the science segment
      % was too short. If it is empty we continue to the next file
      if isempty(gTmp)
        continue;
      end;

    end;

    % Add the results from this file to the existing results in column
    gpsTimes = [ gpsTimes; gTmp];
    naiveSigmas = [ naiveSigmas; nsTmp ];
    ccSigmas = [ ccSigmas; ccsTmp ];
  end;

  % ensure the output is sorted by GPS time
  [ gpsTimes, ind ] = sort(gpsTimes);
  naiveSigmas = naiveSigmas(ind);
  ccSigmas = ccSigmas(ind);

return;
