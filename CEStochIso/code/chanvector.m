function [vector,sampRate,vectorError] = chanvector(chanStr,gpsRange,varargin)
% CHANVECTOR_S5 return vector sub-section of IGWD channel
% [vector[,sampRate,vectorError]]=chanvector(chanStr,gpsStartSec,durationSec)
%    vector is time-series of structure 'chanStr' starting at gpsStartSec
%       with a duration in durationSec.
%
% [vector[,sampRate,vectorError]] = chanvector(chanStr,gpsStartSec:gpsEndSec)
%    vector is time-series of structure 'chanStr' starting at gpsStartSec
%       and ending just before gpsEndSec.
%
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec,durationSec,gpsTimeTbl,frameTbl[,durTbl])
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec:gpsEndSec,gpsTimeTbl,frameTbl[,durTbl])
%     use pre-fetched lists of GPS start times and frame files
%     durTbl = [OPTIONAL] pre-fetched list of frame file durations 
%                        (assume 16 seconds if missing)
%
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec:gpsEndSec,frameListStruct)
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec,durationSec,frameListStruct)
%     use pre-fetched structure of GPS start times,frame files,durations
%
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec:gpsEndSec,lalCacheFile)
% [vector[,sampRate,vectorError]] = chanvector(chanStr,...
%                   gpsStartSec,durationSec,lalCacheFile)
%     gets frame file list from lalCacheFile from LSCdataFind
%
%  *** if inputting frame list information, 'chanStr' can be channel name
%       (i.e. 'H1:LSC-AS_Q') instead of Channel structure
%
% vector would be [] if failure
% sampRate [OPTIONAL] is # of samples/second of frames read in
% vectorError [OPTIONAL] error code from vector creation
%      = 0 if no problems
%        1-10 error from FRAMELIST
%        11-30 error from MLFRAMEQUERY (via FRAMELIST)
%        31-40 error from FETCHSERIES
%        41-50 error within CHANVECTOR
%      41 Channel input is not a structure
%      42 non-numeric GPS range
%      43 input structure is not a Channel structure
%      44 Bad GPS range
%      45 External frame list is empty
%
% $Id: chanvector.m,v 1.4 2009-02-06 14:37:23 jromano Exp $

% IF invalid # of inputs, outputs
%   flag as error
%   exit
% ENDIF
% IF second input is not numeric
%   flag as error
%   exit
% ENDIF
% IF status output requested
%   SET status output = fail (to start)
% ENDIF
error(nargchk(2,6,nargin),'struct');
error(nargchk(1,3,nargout),'struct');
vector = [];
sampRate = 0;
vectorError = 0;
if(~isa(gpsRange,'numeric'))
    vectorError = 42;
    msgId = 'chanvector:badGps';
    error(msgId,'%s: GPS range not numeric',msgId);
end

% - determine if frame list info is expected
% IF Only 2 inputs 
%   CLEAR frame list info flag
% ELSEIF 3 inputs
%   IF third input is numeric
%       CLEAR frame list info flag
%   ELSE
%       SET frame list info flag
% ELSE (4 or more inputs)
%   SET frame list info flag
% ENDIF
frameListFlag = false;
if(nargin <= 2)
    frameListFlag = false;
elseif(nargin ==3)
    inputCand = varargin{1};
    if(isa(inputCand,'numeric'))
        frameListFlag = false;
    else
        frameListFlag = true;
    end
else
    frameListFlag = true;
end

% IF first input is structure
%   get Channel structure properties
%   IF any are empty (i.e. not defined)
%       flag as non-channel struct error
%       EXIT
%   ENDIF
% ELSE 
%   IF no frame list info
%       CREATE error about missing Channel structure
%       EXIT
%   ENDIF
%   IF first input is not a channel string
%       CREATE error about invalid channel input
%       EXIT
%   ENDIF
%   SET channel name from input
% ENDIF
if(isstruct(chanStr))
    if(~isfield(chanStr,'name') || ...
       ~isfield(chanStr,'instrument') || ...
       ~isfield(chanStr,'type') || ...
       ~isfield(chanStr,'site') || ...
       ~isfield(chanStr,'rate'))
        vectorError = 43;
        msgId = 'chanvector:badChannel';
        error(msgId,'%s: Channel input is not a Channel structure',msgId);
    end
    name = chanStr.name;
    instrument = chanStr.instrument;
    type = chanStr.type;
    site = chanStr.site;
    rate = chanStr.rate;
    if(isempty(name) || isempty(type) || isempty(site) ||...
         isempty(instrument) || isempty(rate))
        vectorError = 43;
        msgId = 'chanvector:badChannel';
        error(msgId,'%s: Channel input is not a Channel structure',msgId);
    end
else
    if(frameListFlag == false)
        vectorError = 41
        msgId = 'chanvector:badChannel';
        error(msgId,'%s: Channel input must be structure',msgId)
    else
        if(~ischar(chanStr))
            vectorError = 41
            msgId = 'chanvector:badChannel';
            error(msgId,'%s: Channel input is not structure or string',msgId);
        end
        colonPos = strfind(chanStr,':');
        if(isempty(colonPos))
            vectorError = 43
            msgId = 'chanvector:badChannel';
            error(msgId,'%s: Channel name %s not valid',msgId,chanStr);
        end
        name = chanStr;
        rate = 0;
        instrument = chanStr(1:1);
        site = [];
        type = [];
    end
end
        
% CLEAR GPS range flag
% CLEAR frame input flag
% CLEAR duration input flag
% IF # columns in GPS Range entry = 1 (indicating #,#)
%   IF there is a duration entry
%       IF it is a valid number
%           SET gpsStartSec to first entry
%           SET durationSec to second entry
%           SET GPS range flag
%           IF pre-fetched GPS start, frames were input
%               SET GPS start time table from input
%               SET frame table from input
%               SET frame input flag
%               IF pre-fetched frame durations were input
%                   SET duration table from input
%                   SET duration input flag
%               ENDIF
%           ENDIF
%       ELSE
%           CREATE error about bad duration value
%           EXIT
%       ENDIF
%   ELSE
%       CREATE error about missing duration
%       EXIT
%   ENDIF
isGpsRange = false;
framesInput = false;
durInput = false;
[mrows,ncols] = size(gpsRange);
if(ncols == 1)
    if(nargin > 2)
        if(isa(varargin{1},'numeric'))
            gpsOffsetSec = gpsRange;
            durationSec = varargin{1};
            isGpsRange = true;
            if(frameListFlag == true)
                if(isa(varargin{2},'numeric'))
                    if(nargin > 4)
                        framesInput = true;
                        gpsTimeTbl = varargin{2};
                        frameTbl = varargin{3};
                        if(nargin > 5)
                            durTbl = varargin{4};
                            durInput = true;
                        end
                    else
                        msgId = 'chanvector:badList';
                        error(msgId,...
                            '%s: Need GPS times and frame list input',msgId);
                    end
                else
                    framesInput = true;
                    listTry = varargin{2};
                    [listError,gpsTimeTbl,frameTbl,durTbl] = ...
                        getframelistinfo(listTry);
                    if(listError ~= 0)
                        msgId = 'chanvector:badframelist';
                        warning(msgId,'%s: frame list input error',...
                            msgId,listError);
                    end
                    if(~isempty(durTbl))
                        durInput = true;
                    end
                end
            else
                framesInput = false;
            end
        else
            vectorError = 44;
            msgId = 'chanvector:badDuration';
            error(msgId,'%s: duration input is not numeric',msgId);
        end
    else
        vectorError = 44;
        msgId = 'chanvector:badDuration';
        error(msgId,'%s: range input has GPS Start but no duration',msgId);
    end

% ELSE
%   SET gpsOffsetSec to 1st element of input vector
%   SET gpsEndSec to last element of input vector
%   SET durationSec = gpsEndSec - gpsOffsetSec
%   SET GPS range flag
%   IF pre-fetched GPS start, frames were input
%       SET GPS start time table from input
%       SET frame table from input
%       SET frame input flag
%       IF pre-fetched frame durations were input
%           SET duration table from input
%           SET duration input flag
%       ENDIF
%   ENDIF
% ENDIF
else
    gpsOffsetSec = gpsRange(1);
    gpsEndSec = gpsRange(ncols);
    durationSec = gpsEndSec - gpsOffsetSec;
    isGpsRange = true;
    if(frameListFlag == true)
        if(isa(varargin{1},'numeric'))
            if(nargin > 3)
                framesInput = true;
                gpsTimeTbl = varargin{1};
                frameTbl = varargin{2};
                if(nargin > 4)
                    durTbl = varargin{3};
                    durInput = true;
                end
            else
                msgId = 'chanvector:badList';
                error(msgId,...
                       '%s: Need GPS times and frame list input',msgId);
            end
        else
            framesInput = true;
            listTry = varargin{1};
            [listError,gpsTimeTbl,frameTbl,durTbl] = ...
                        getframelistinfo(listTry);
            if(listError ~= 0)
                msgId = 'chanvector:badframelist';
                warning(msgId,'%s: frame list input error',msgId,listError);
            end
            if(~isempty(durTbl))
                durInput = true;
            end
        end
    else
       framesInput = false;
    end
end

% IF gps Range flag not set
%   FLAG as error
%   exit
% ENDIF
% IF duration is not > 0
%   FLAG as error
%   exit
% END
% IF channel input as object
%   SET start time = channel start time + GPS offset
% ELSE
%   SET start time = GPS offset
% ENDIF
if(isGpsRange ~= true)
    vectorError = 44;
    error('CHANVECTOR: Bad GPS range format');
end
if(durationSec <= 0)
    vectorError = 44;
    error('CHANVECTOR: GPS range duration is <= 0');
end
if(isstruct(chanStr))
    startTime = chanStr.gpsStart + gpsOffsetSec;
else
    startTime = gpsOffsetSec;
end

% IF frame input flag is cleared (the default)
%   USE framelist function to get list of frame files covering
%            start time to (start time + duration)
%   IF framelist was not successful
%       FLAG as error
%   ENDIF
%   SET duration input flag
% ELSE
%   IF external frame lists are empty
%       CREATE error
%       EXIT
%   ELSE
%       IF length of GPS,frame,durations do not match
%           CREATE error
%           exit
%       ELSE
%           SORT the input frame file list by starting GPS
%           REMOVE duplicate files in lists
%       ENDIF
% ENDIF
if(framesInput ~= true)
    [gpsTimeTbl,frameTbl,listError,durTbl] = framelist(instrument,...
            type,startTime,durationSec);
    vectorError = listError;
    if(listError ~= 0)
        msgId = 'chanvector:badframelist';
        errStr = ...
            sprintf('CHANVECTOR: error %d in getting frame file list\n',...
            listError);
        warning(msgId,errStr);
    end
    durInput = true;
else
    nTimes = length(gpsTimeTbl);
    nFiles = length(frameTbl);
    if(true == durInput)
        nDurs = length(durTbl);
    else
        nDurs = nTimes;
    end
    vectorError = 0;
    if(nTimes < 1 || nFiles < 1 || nDurs < 1)
        tmp = sprintf('CHANVECTOR: external frame lists are empty\n');
        vectorError = 45;
        error(tmp);
    else
        if((nTimes ~= nFiles) || (nDurs ~= nFiles))
            tmp = sprintf(...
      'CHANVECTOR: external frame,GPS,duration lists do not match!\n');
            vectorError = 46;
            error(tmp);
        else
            [gpsTimeTbl,uniqueIndx]=unique(gpsTimeTbl);
            frameTbl = {frameTbl{uniqueIndx}};
            if(true == durInput)
                durTbl = durTbl(uniqueIndx);
            end
        end
    end
end

% IF duration list exists
%   Use fetchseries with duration list input to create series from file list
% ELSE
%   Use fetchseries without duration list input
% END
% SET channel rate to frame rate
% IF series has missing data
%   CREATE error string
% ENDIF
% IF vector status requested
%   SET output to fetchseries output flag
% ENDIF
if(durInput == true)
    [vector,trueRate,seriesError] = fetchseries(name,rate,startTime,...
        durationSec,gpsTimeTbl,frameTbl,durTbl);
else
    [vector,trueRate,seriesError] = fetchseries(name,rate,startTime,...
        durationSec,gpsTimeTbl,frameTbl);
end
sampRate = trueRate;
if(seriesError ~= 0)
    msgId = 'chanvector:badtimeseries';
    if (0 == vectorError)
        vectorError = seriesError + 30;
    end
    tmp = sprintf('CHANVECTOR: Error %d when retrieving time-series!\n',...
        seriesError);
    warning(msgId,tmp);
end
return
