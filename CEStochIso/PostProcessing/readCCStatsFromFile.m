function [gpsTimes, ccStats, ccSigmas] = readCCStatsFromFile(fileList, isComplex)
%
%  readCCStatsFromFile -- read naive sigmas from file or list of files
%
%  [ gpsTimes, ccStats, ccSigmas ] = readCCStatsFromFile(fileList, isComplex)
%
%  returns column vectors containing GPS times, CC stats and CC sigmas
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
%    <gpsTime> <CC> <ccSigma>
%
%  or Matlab files produced by stochastic.m using the writeOutputToMatFile flag.
%
%  The flag isComplex is required only for text files and tells the function
%  that the CC statistic is complex and consists of two number, the real
%  and the imaginary part. If it is omitted the CC statistic is assumed to
%  be real.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (nargin < 2)
    isComplex = false;
  end;

  gpsTimes = [];
  ccStats = [];
  ccSigmas = [];

  % If the file list is just a single filename, convert
  % it to a cellarray with one element.
  if (~iscell(fileList))
    fileList = { fileList };
  end;

  for f = fileList
    f = f{1};
    checkFileExists('CC stats file', f);

    [pathstr, name, ext] = fileparts(f);
    if (ext == '.mat');
      % Check if there are no intervals to include, which can happen if eg. the science segment
      % was too short. If there are no segments we continue to the next file
      load(f, 'params');
      if (params.numIntervalsTotal < 1)
        continue;
      end;

      % Matfile doesn't care if ccStat is real or complex so we can just load it
      load(f, 'segmentStartTime', 'ccStat', 'ccSigma');
      gTmp = segmentStartTime;
      % Stats are stored in the .mat file as ccSigma(<segment index>, <trial index>). Post processing
      % doesn't handle multiple trials so we always use trial 1.
      ccTmp = ccStat(:, 1);
      ccsTmp = ccSigma(:, 1);
      clear segmentStartTime ccStat ccSigma;
    else
      % Read from text file
      if (isComplex)
        [gTmp, ccTmpReal, ccTmpImag, ccsTmp ] = textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');
        ccTmp = complex(ccTmpReal, ccTmpImag);
      else
        [gTmp, ccTmp, ccsTmp] = textread(f, '%f%f%f\n', -1, 'commentstyle', 'matlab');
      end;

      % Check for empty results, which can happen if eg. the science segment
      % was too short. If it is empty we continue to the next file
      if isempty(gTmp)
        continue;
      end;

    end;

    gpsTimes = [ gpsTimes; gTmp];
    ccStats = [ ccStats; ccTmp ];
    ccSigmas = [ ccSigmas; ccsTmp ];
  end;

  % ensure the output is sorted by GPS time
  [gpsTimes, ind] = sort(gpsTimes);
  ccStats = ccStats(ind);
  ccSigmas = ccSigmas(ind);

return;
