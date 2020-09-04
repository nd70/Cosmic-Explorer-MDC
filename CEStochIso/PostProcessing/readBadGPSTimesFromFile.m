function gpsTimes = readBadGPSTimesFromFile(fileList)
%
%  readBadGPSTimesFromFile -- read the bad GPS times from a file or list of files
%
%  gpsTimes = readBadGPSTimesFromFile(fileList)
%
%  returns column vector containing all the bad GPS times
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
%    <gpsTime>
%
%  or Matlab files produced by stochastic.m using the writeOutputToMatFile flag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gpsTimes = [];

  % If the file list is just a single filename, convert
  % it to a cellarray with one element.
  if (~iscell(fileList))
    fileList = { fileList };
  end;

  for f = fileList
    f = f{1};
    checkFileExists('Bad GPS times file', f);

    [pathstr, name, ext] = fileparts(f);
    if (ext == '.mat');

      % Check if there are no intervals to include, which can happen if eg. the science segment
      % was too short. If there are no segments we continue to the next file
      load(f, 'params');
      if (params.numIntervalsTotal < 1)
        continue;
      end;

      load(f, 'badGPSTimes');
      gTmp = badGPSTimes{1};
    else
      % Read from text file
      gTmp = textread(f, '%f\n', -1, 'commentstyle', 'matlab');

      % Check for empty results, which can happen if eg. the science segment
      % was too short. If it is empty we continue to the next file
      if isempty(gTmp)
        continue;
      end;

    end;

    % Add the results from this file to the existing results in column
    gpsTimes = [ gpsTimes; gTmp];
  end;

  % ensure the output is sorted by GPS time
  [ gpsTimes, ind ] = sort(gpsTimes);

return;
