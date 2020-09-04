function [outData,rate,seriesError] = ...
    fetchseries(channel,defRate,gpsStart,duration,...
		gpsTimes,frameFiles,frameDurs)
% FETCHSERIES - builds time-series for channel from GPS interval,
% list of frame files
% 
% loads channel data from frame files using frgetvect
% and returns the time-series data, sample rate, status flag
%
% [vector,rate,seriesError] = ...
%    fetchseries(channel,defRate,gpsStart,duration,...
%                gpsTimes,frameFiles[,frameDurs]);
% channel is the name of the desired channel and should be given as a string
%       example: channel = 'H2:LSC-AS_Q'
% defRate is a default sample rate (can be 0)
% gpsStart is the GPS start time (whole seconds).  It does not have
%       to correspond to the start time of any frame file.
% duration is the length of the desired data segment. Its value
%       should be specified in whole seconds
% gpsTimes - list of starting GPS time for each frame file
% frameFiles - list of paths for frame files to be read in order
% frameDurs  - (optional) list of durations for each frame file
%   NOTE: frame files MUST be sorted in ascending GPS order!
%
%    vector = time-series
%    rate = [OPTIONAL] sample rate of time series
%    seriesError = [OPTIONAL] error code = 0 if no problem, >0 otherwise
%       1 no files in list
%       2 list doesn't overlap time range
%       3 file in list not found
%       4 files do not completely cover time range
%       5 time-series data in files has 0's
%       6 unable to determine rate
%       7 no data for channel
%
% Warnings:
%   fetchseries:emptyFrameList
%   fetchseries:missingFiles
%   fetchseries:rateUndetermined
%   fetchseries:fileNotFound
%   fetchseries:frgetvect
%   fetchseries:zeroData
%   fetchseries:missingData
%   fetchseries:rateChanged
%   fetchseries:noData
%
% Errors:
%   fetchseries:inOutArg
%   fetchseries:emptyFrameList
%   fetchseries:missingFiles
%   fetchseries:rateUndetermined
%   fetchseries:fileNotFound
%   fetchseries:frgetvect
%   fetchseries:zeroData
%   fetchseries:missingData
%   fetchseries:noData
%
% $Id: fetchseries.m,v 1.7 2009-02-06 14:37:23 jromano Exp $

% IF the # of arguments is wrong
%    SET failed status
%    CREATE error print showing calling syntax
%    RETURN
% ENDIF
inErr = nargchk(6,7,nargin);
outErr = nargoutchk(1,3,nargout);
if(~isempty(inErr) || ~isempty(outErr))
    msgId = 'fetchseries:inOutArg';
    errmsg = sprintf(...
      '%s:\n\t%s...\n\t%s',...
      msgId,...
      '[data[,rate,OK]] = fetchseries(channel,defRate,gpsStart,',...
      'duration,gpsTimes,frameFiles[,frameDurs])');
    error(msgId,errmsg);
end
errFlag = nargout < 3;
DEFAULT_FRAME_SIZE = 16;

% - added to turn off informative warning from frgetvect
warning('off','frgetvect:info');

% CALCULATE series end time
% CALCULATE GPS start second, GPS end second
% IF sample rate is pre-defined for the channel
%   CALCULATE true start and duration to ensure whole samples
%   SET output series = (series duration * sample rate) samples = 0
% ELSE
%   CLEAR sample rate
%   CLEAR output series (i.e. make it a null)
% ENDIF
outData = [];
rate = 0;
seriesError = 0;
gpsEnd = gpsStart + duration;
[startWholeSec,endWholeSec,trueStart,trueDur] = ...
    calcseriesrange(gpsStart,duration,defRate);
if(defRate > 0)
    rate = defRate;
    if(trueStart > 0 && trueDur > 0)
        trueEnd = trueStart + trueDur;
        outData = zeros([1,trueDur*rate]);
    else
        trueEnd = 0;
        outData = [];
    end
else
    rate = 0;
    trueSec = 0;
    trueDur = 0;
    trueEnd = 0;
    outData = [];
end

% IF no frame files were input
%   OUTPUT error or warning for no input files
%   SET output series to empty set
%   SET sampling rate = 0
%   EXIT
% ENDIF
% CALCULATE Frame End Times
% IF frame file durations were input
%   SET file lengths to input
% ELSE
%   SET file lengths to default
% ENDIF
numFiles = length(gpsTimes);
if(numFiles < 1)
  msgId = 'fetchseries:emptyFrameList';
  if(errFlag)
    error(msgId,msgId);
  else
    seriesError = 1;
    warning(msgId,msgId);
    return
  end
end
if(nargin > 6)
    frameSecs = frameDurs;
else
    frameSecs(1:numFiles) = DEFAULT_FRAME_SIZE;
end
endTimes = gpsTimes + frameSecs;

% FIND frames with end times after GPS start second
% IF such frames were found
%   SET index of first used frame to first frame in list
% ELSE
%   OUTPUT error about lack of files
%   EXIT
% ENDIF
findx = find(endTimes > startWholeSec);
if(length(findx) > 0)
  frstFile = findx(1);
else
    seriesError = 2;
    msgId = 'fetchseries:missingFiles';
    warning(msgId,'%s: requested interval [%d, %d]',...
	  msgId,gpsStart,gpsEnd);
    for ii=1:numFiles
      warning(msgId,...
	      '%s:frame file %d has GPS start %d duration %d end %d',...
	      msgId,ii,gpsTimes(ii),frameSecs(ii),endTimes(ii));
    end
    if(errFlag)
        error(msgId,msgId);
    end
    return
end

% FIND frames with start times before GPS end second
% IF such frames were found
%   SET index of last used frame to last frame in list
% ELSE
%   OUTPUT error about lack of files
%   EXIT
% ENDIF
lstIndx = find(gpsTimes < endWholeSec);
if(length(lstIndx) > 0)
    lastFile = lstIndx(end);
else
    seriesError = 2;
    msgId = 'fetchseries:missingFiles';
    warning(msgId,'%s: requested interval [%d, %d]',...
	  msgId,gpsStart,gpsEnd);
    for ii=1:numFiles
        warning(...
	      msgId,'%s: frame %d has GPS start %d duration %d end %d \n',...
	      msgId,ii,gpsTimes(ii),frameSecs(ii),endTimes(ii));
    end
    if(errFlag)
        error(msgId,msgId);
    end
    return
end

%  -- retrieve needed data from the first frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF GEO FRAME >>>>>>>>>>>>
%keyboard
if strncmp(channel,'G',1)
  % && May, 2008 - To get around problems with reading partial frames,
  %   we will read the whole frame for the channel, then
  %  hack out the piece we need &&
  % SET read limits for whole frame
  % IF first frame file does not exist
  %   CREATE Error about missing file
  %   EXIT
  % ENDIF
  % GET data from first frame
  frstStart = gpsTimes(frstFile);
  frstEnd = endTimes(frstFile);
else
  % SET starting time = Maximum (first frame start time, series start second)
  % SET ending time = Minimum (first frame end time, series end second)
  % SET duration in frame = ending time - starting time
  % IF first frame file does not exist
  %   CREATE Error about missing file
  %   EXIT
  % ENDIF
  % SET temp vector = data from first frame from starting time to ending time
  % CALCULATE sampling rate = vector length / frame duration
  frstStart = max(gpsTimes(frstFile),startWholeSec);
  frstEnd = min(endTimes(frstFile),endWholeSec);
end
%% IF GEO FRAME <<<<<<<<<<<
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frstDur = frstEnd - frstStart;
fileStr = char(frameFiles{frstFile});
if(exist(fileStr,'file') == 0)
    fileNotFound(fileStr,errFlag);
    seriesError = 3;
    return
end
%keyboard
try
    tmpVect = frgetvect(fileStr,channel,frstStart,frstDur);
catch
    tmpVect = [];
    errmsg = lasterr;
    frgetvectErr(errmsg,fileStr,channel,frstStart,frstDur,errFlag);
end
% CALCULATE sampling rate = vector length / frame duration
%  trap cases where no data is returned 
%       (typically for a non-existent channel)
if(length(tmpVect)<1)
    msgId = 'fetchseries:noData';
    seriesError = 7;
    if(errFlag)
        error(msgId,'%s: Channel %s has no data from frame %s',...
            msgId,channel,fileStr);
    else
        warning(msgId,'%s: Channel %s has no data from frame %s',...
            msgId,channel,fileStr);
    end    
    rate = 0;
else
    rate = length(tmpVect)/frstDur;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF GEO FRAME >>>>>>>>>>>>
if strncmp(channel,'G',1)

  % - hack out the part of the data vector we need
  % IF starting time > start of first frame
  %   CALCULATE first sample = 1 + time-diff/(rate in Hz)
  % ELSE
  %   SET first sample = 1
  % ENDIF
  % IF ending time < end of first frame
  %   CALCULATE last sample = first sample + duration/(rate in Hz)
  % ELSE
  %   SET last sample = last sample in frame
  % ENDIF
  % SET temp vector = input vector from updated first to last sample
  begSamp = 1;
  endSamp = length(tmpVect);
  startDiff = startWholeSec - frstStart;
  if(startDiff > 0)
      frstSamp = round(startDiff*rate) + begSamp;
  else
      frstSamp = begSamp;
  end
  endDiff = frstEnd - endWholeSec;
  if(endDiff > 0)
      timeDiff = endWholeSec - startWholeSec;
      sampDiff = round(timeDiff*rate);
      lastSamp = frstSamp + sampDiff;
  else
      lastSamp = endSamp;
  end
  tmpVect = tmpVect(frstSamp:lastSamp);
  frstStart = max(gpsTimes(frstFile),startWholeSec);
  frstEnd = min(endTimes(frstFile),endWholeSec);

end
%% IF GEO FRAME <<<<<<<<<<<
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -- Here we re-do the true start, end time calculation
%  as we have an 'official' rate from a real frame file.
%
% IF rate is still 0
%   IF pre-set rate is non-zero
%       SET rate from pre-set rate
%   ELSE
%       CREATE error that rate is unknown
%       EXIT
%   ENDIF
% ELSE
%   IF the pre-set rate was non-zero
%       IF true rate not the same as preset rate
%           CREATE warning
%       ENDIF
%   ENDIF
% ENDIF
% CALCULATE true starting time, duration to align with samples
% CALCULATE true ending time = true starting time + true duration
% INITIALIZE output data-series to correct length, but all zeros
if (rate <= 0)
    if defRate > 0
        rate = defRate;
    else
        msgId = 'fetchseries:rateUndetermined';
        if(seriesError < 1)
            seriesError = 6;
        end
        if(errFlag)
            error(msgId,msgId);
        else
            warning(msgId,msgId);
        end
    end
else
    if defRate > 0
        if(defRate ~= rate)
            msgId = 'fetchseries:rateChanged';
            msg = sprintf('True rate %f does not match pre-set rate %f\n',rate,defRate);
            warning(msgId,msg);
        end
    end
end
[startWholeSec,endWholeSec,trueStart,trueDur] = ...
    calcseriesrange(gpsStart,duration,rate);
trueEnd = trueStart + trueDur;
outData = zeros([1,trueDur*rate]);

% SET first sample = ...
%   (true start - starting time used for first frame) * sample rate
% IF 'first sample' < 1 (series start is before data retrieved)
%   OUTPUT warning about missing data
%   SET bad data flag
%   SET next output sample = 2 - 'first sample'
%   SET 'first sample' = 1
% ELSE
%   SET next output sample = 1
% ENDIF
frstSamp = (trueStart - frstStart)*rate + 1;
%D fprintf('FETCHSERIES: original first sample is %g\n',frstSamp);
% clamp to integer
frstSamp = round(frstSamp);
if(frstSamp < 1)
    if(seriesError < 1)
        seriesError = 4;
    end
    missingData(channel,trueStart,abs(frstSamp)/rate+trueStart,errFlag);
    nextOut = 2 - frstSamp;
    frstSamp = 1;
else
    nextOut = 1;
end

% IF first frame not same as last frame (at least two frames for series)
%   IF data from first frame is non-zero
%       SET output series from 'next output' on = temp vector from 'first sample' to end of vector
%   ELSE
%       CREATE warning error about data series = 0
%       SET error flag
%   ENDIF
%   IF first frame index - last frame index > 1 (there are full frames in the middle)
%       LOOP from first index+1 to last index-1
%           IF frame file does not exist
%               CREATE error message
%               EXIT
%           ENDIF
%           SET temp vector = time-series from frame file
%           IF data is all zeros
%               OUTPUT warning about bad data
%               SET bad data flag
%           ELSE
%               SET gap = frame's start time - (previous frame's start time + previous frame duration)
%               IF gap > 0 (indicating a gap in data)
%                   OUTPUT warning about missing data
%                   SET bad data flag
%               ENDIF
%               SET next output sample = (frame start time - series start)*sample rate + 1
%               SET output series from 'next output' to 'next output' + vector length = temp
%           ENDIF
%       ENDLOOP
%   ENDIF
if(frstFile < lastFile)
    sampSpan = length(tmpVect)-frstSamp;
    lastOut = nextOut+sampSpan;
    if(iszero(tmpVect))
        seriesError = 5;
        zeroData(channel,fileStr,frstStart,frstEnd,errFlag);
    else
        outData(nextOut:lastOut)=tmpVect(frstSamp:end);
    end
    if((lastFile-frstFile)>1)
        for iFile=frstFile+1:lastFile-1
            fileStr = char(frameFiles{iFile});
            if(exist(fileStr,'file') == 0)
                fileNotFound(fileStr,errFlag);
                seriesError = 3;
                return
            end
            try
                tmpVect = ...
                  frgetvect(fileStr,channel,...
                  gpsTimes(iFile),frameSecs(iFile));
            catch
                tmpVect = [];
                errmsg = lasterr;
                frgetvectErr(errmsg,fileStr,channel,...
                  gpsTimes(iFile),frameSecs(iFile),errFlag);
	    end
	    if(length(tmpVect)<1)
            msgId = 'fetchseries:noData';
            seriesError = 7;
            if(errFlag)
                error(msgId,'%s: Channel %s has no data from frame %s',...
                    msgId,channel,fileStr);
            else
                warning(msgId,'%s: Channel %s has no data from frame %s',...
                    msgId,channel,fileStr);
            end
        end
        if(iszero(tmpVect))
	      if(seriesError < 1)
              seriesError = 5;
          end
          zeroData(channel,fileStr,gpsTimes(iFile),...
		       gpsTimes(iFile)+frameSecs(iFile),errFlag);
	    else
	      gapSec = gpsTimes(iFile) - ...
		       (gpsTimes(iFile-1) + frameSecs(iFile-1));
	      if(gapSec > 0)
            if(seriesError < 1)
                seriesError = 4;
            end
            missingData(channel,...
		            (gpsTimes(iFile-1)+frameSecs(iFile-1)),...
		            gpsTimes(iFile),errFlag);
	      end
	      nextOut = (gpsTimes(iFile)-trueStart)*rate + 1;
%D	      fprintf('FETCHSERIES: original nextOut is %g\n',nextOut);
% clamp index to integer value
          nextOut = round(nextOut);
          lastOut = nextOut+length(tmpVect)-1;
	      outData(nextOut:lastOut) = tmpVect;
	    end
        end
    end

    %  --- get data from last frame file (when it isn't also the first frame file)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IF GEO FRAME >>>>>>>>>
    if strncmp(channel,'G',1)
      % && May, 2008 - To get around problems with reading partial frames,
      %   we will read the whole frame for the channel, then
      %  hack out the piece we need &&
      %
      %   IF frame file does not exist
      %       CREATE error message
      %       EXIT
      %   ENDIF
      %   SET temp vector = time series from last frame used
      lastStart = gpsTimes(lastFile);
      lastEnd = endTimes(lastFile);
    else
      %
      %   IF frame file does not exist
      %       CREATE error message
      %       EXIT
      %   ENDIF
      %   SET temp vector = time series from last frame used
      %   IF data is all zeros
      %       OUTPUT warning for bad data
      %       SET bad data flag
      %   ELSE
      %       SET gap = last frame's start time - (previous frame's start time + previous frame duration)
      %       IF gap > 0 (indicating a gap in data)
      %           OUTPUT warning about missing data
      %           SET bad data flag
      %       ENDIF
      %       SET first sample = beginning of vector
      %       SET next output sample = (frame start time - series start)*sample rate + 1
      %   ENDIF
      % ELSE (only one frame file)
      %   SET starting second in last frame = that in first frame
      %   SET ending second in last frame = that in first frame
      %   SET duration in last frame = that in first frame
      % END
      lastStart = gpsTimes(lastFile);
      lastEnd = min(endTimes(lastFile),endWholeSec);
    end
    % IF GEO FRAME <<<<<<<<<<<<<
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lastDur = lastEnd - lastStart;
    fileStr = char(frameFiles{lastFile});
    if(exist(fileStr,'file') == 0)
        fileNotFound(fileStr,errFlag);
        seriesError = 3;
        return
    end
    try
        tmpVect = frgetvect(fileStr,channel,lastStart,lastDur);
    catch
        tmpVect = [];
        errmsg = lasterr;
        frgetvectErr(errmsg,fileStr,channel,lastStart,lastDur,errFlag);
    end
    if(length(tmpVect)<1)
        msgId = 'fetchseries:noData';
        seriesError = 7;
        if(errFlag)
            error(msgId,'%s: Channel %s has no data from frame %s',...
                msgId,channel,fileStr);
        else
            warning(msgId,'%s: Channel %s has no data from frame %s',...
                msgId,channel,fileStr);
        end    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % IF GEO FRAME >>>>>>>>>
    if strncmp(channel,'G',1)
      % - hack out the part of the data vector we need
      %   SET first sample = 1
      %   IF ending time < end of first frame
      %       CALCULATE last sample = last sample in frame - time-diff/rate
      %   ELSE
      %       SET last sample = last sample in frame
      %   ENDIF
      %   SET temp vector = input vector from updated first to last sample
      begSamp = 1;
      endSamp = length(tmpVect);
      frstSamp = begSamp;
      endDiff = lastEnd - endWholeSec;
      if(endDiff > 0)
          lastSamp = endSamp - round(endDiff*rate);
      else
          lastSamp = endSamp;
      end
      tmpVect = tmpVect(frstSamp:lastSamp);
      lastEnd = min(endTimes(lastFile),endWholeSec);
    
      %   IF data is all zeros
      %       OUTPUT warning for bad data
      %       SET bad data flag
      %   ELSE
      %       SET gap = last frame's start time - (previous frame's start time + previous frame duration)
      %       IF gap > 0 (indicating a gap in data)
      %           OUTPUT warning about missing data
      %           SET bad data flag
      %       ENDIF
      %       SET first sample = beginning of vector
      %       SET next output sample = (frame start time - series start)*sample rate + 1
      %   ENDIF
      % ELSE (only one frame file)
      %   SET starting second in last frame = that in first frame
      %   SET ending second in last frame = that in first frame
      %   SET duration in last frame = that in first frame
      % END
    end
    % IF GEO FRAME <<<<<<<<<<<<<
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(iszero(tmpVect))
        if(seriesError < 1)
            seriesError = 5;
        end
        zeroData(channel,fileStr,lastStart,lastEnd,errFlag);
    else
        gapSec = gpsTimes(lastFile) - (gpsTimes(lastFile-1) ...
                                     + frameSecs(lastFile-1));
        if(gapSec > 0)
                if(seriesError < 1)
                seriesError = 4;
            end
            missingData(channel,...
                  (gpsTimes(lastFile-1)+frameSecs(lastFile-1)),...
                  gpsTimes(lastFile),errFlag);
        end
        frstSamp = 1;
        nextOut = (gpsTimes(lastFile)-trueStart)*rate + 1;
        %D          fprintf('FETCHSERIES: original nextOut is %g\n',nextOut);
        %  clamp index to integer value
        nextOut = round(nextOut);
    end
else
  lastStart = frstStart;
  lastEnd = frstEnd;
  lastDur = frstDur;
end

% SET last sample = (series end - starting time in last frame) * sample rate
% IF 'last sample' exceeds temp vector length (series end is after frame)
%   OUTPUT warning about missing data
%   SET bad data flag
%   SET last sample = end of vector
% END
% SET output series = output series + temp vector from 'first sample' to 'last sample'
% SET time series output = accumulated output series
% IF error flag output requested
%   SET error flag output
% ENDF
lastSamp = (trueEnd - lastStart) * rate;
%D fprintf('FETCHSERIES: original last sample is %g\n',lastSamp);
% clamp index to integer value
lastSamp = round(lastSamp);
if(lastSamp > length(tmpVect))
    if(seriesError < 1)
        seriesError = 4;
    end
    missingData(channel,gpsTimes(lastFile),trueEnd,errFlag);
    lastSamp = length(tmpVect);
end
sampSpan = lastSamp - frstSamp;
lastOut = nextOut+sampSpan;
%D fprintf('FETCHSERIES: frstSamp %g lastSamp %g sampSpan %g nextOut %g lastOut %g\n',...
%D    frstSamp,lastSamp,sampSpan,nextOut,lastOut);
if(~iszero(tmpVect))
    outData(nextOut:lastOut)=tmpVect(frstSamp:lastSamp);
end
return

function result = iszero(dataSeries)
% ISZERO - tests for vector with all elements = 0;
%
%   iszero(dataSeries)
%       dataSeries = any vector, matrix
%    iszero = 1 if empty or all values = 0
%    iszero = 0 otherwis

if(~isempty(dataSeries))
    minVal = min(dataSeries);
    maxVal = max(dataSeries);
    if (minVal == 0 ) && (maxVal == 0)
        result = 1;
    else
        result = 0;
    end
else
    result = 1;
end
return

function fileNotFound(fileStr,errFlag)
% FILENOTFOUND - signal error or warning if file not found
%
% fileNotFound(fileStr,errFlag)
%
% fileStr - name of file
% errFlag - true if error, false if warning

msgId = 'fetchseries:fileNotFound';
if (errFlag)
    error(msgId,'%s: %s',msgId,fileStr);
else
    warning(msgId,'%s: %s',msgId,fileStr);
end
return

function frgetvectErr(errmsg,fileStr,chan,start,dur,errFlag)
% FRGETVECTERR - warning if frgetvect errors out
%
% frgetvectErr(errmsg,fileStr,chan,start,dur,errFlag)
%
% errmsg   - frgetvect error message
% fileStr  - file name where error occured
% chan     - channel being read
% start    - gps start time
% dur      - duration
% errFlag  - true if error, false if warning

msgId = 'fetchseries:frgetvect';
msg = sprintf(...
    '%s:%s file %s, channel %s from %d for %d',...
    msgId,errmsg,fileStr,chan,start,dur);
if(errFlag)
    error(msgId,msg);
else
    warning(msgId,msg);
end
return

function zeroData(chan,fname,gpsStart,gpsEnd,errFlag)
% ZERODATA - signal warning if channel data == 0
%
% zeroData(chan,fname,gpsStart,gpsEnd,errFlag)
%
% chan       - channel name
% fname      - frame file name
% gpsStart   - gps time where 0's start
% gpsEnd     - gps time where 0's end
% errFlag  - true if error, false if warning

msgId = 'fetchseries:zeroData';
msg = sprintf(...
    '%s: %s in %s = 0 over range [%d,%d]',...
    msgId,chan,fname,gpsStart,gpsEnd);
if(errFlag)
    error(msgId,msg);
else
    warning(msgId,msg);
end
return

function missingData(chan,gpsStart,gpsEnd,errFlag)
% MISSINGDATA - warning if data is missing
%
% missingData(chan,start,end,errFlag)
%
% chan       - channel name
% gpsStart   - gps start time of gap
% gpsEnd     - gps end time of gap
% errFlag  - true if error, false if warning

msgId = 'fetchseries:missingData';
msg = sprintf('%s: %s from %d for %d',...
	msgId,chan,...
	gpsStart,gpsEnd-gpsStart);
if(errFlag)
    error(msgId,msg);
else
    warning(msgId,msg);
end
return

