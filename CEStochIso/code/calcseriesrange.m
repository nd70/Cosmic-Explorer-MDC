function [startWholeSec,endWholeSec,trueStart,trueDur] = calcseriesrange(gpsStart,duration,rate)
% CALCSERIESRANGE - calculate time-series limits for FETCHSERIES
% 
% [startWholeSec,endWholeSec,trueStart,trueDur] = calcseriesrange(gpsStart,duration,rate)
%      gpsStart - start time of series (in GPS seconds) can be non-integer
%      duration - duration of time series (in seconds)
%      rate - sample rate of time series (can be 0)
%
%      startWholeSec - starting GPS in whole seconds (for frame file checks)
%      endWholeSec - ending GPS in whole seconds (for frame file checks)
%      trueStart - true starting GPS second ( = 0 if rate = 0)
%      trueDur - true series duration in seconds ( = 0 if rate = 0)
%      
%    true start, duration may differ from input values to match granularity
%       of the sample rate.  This is to guarantee whole # of samples

% $Id: calcseriesrange.m,v 1.1 2004-08-11 18:43:34 kathorne Exp $

% IF the # of arguments is wrong
%    SET failed status
%    CREATE error print showing calling syntax
%    RETURN
% ENDIF
% CALCULATE GPS start second from series start time
% CALCULATE series end time
% CALCULATE GPS end second from series end time
inErr = nargchk(3,3,nargin);
outErr = nargoutchk(4,4,nargout);
if(~isempty(inErr) | ~isempty(outErr))
    fmtStr = '[startWholeSec,endWholeSec,trueStart,trueDur] = calcseriesrange(gpsStart,duration,rate)';
    fprintf(' ERROR - calcseriesrange syntax is %s\n',fmtStr);
    return
end
startWholeSec = floor(gpsStart);
gpsEnd = gpsStart + duration;
endWholeSec = ceil(gpsEnd); 

% IF rate is defined (> 0)
%   IF gps start not in whole seconds
%       GET fractional seconds
%       SET true fraction = [rounded down (fraction *sample rate)] / sample rate
%       SET true start = gps start (whole seconds) + true fraction
%   ELSE
%       SET true start = gps start (whole seconds)
%   ENDIF
%   SET true duration = [rounded down(duration*sample rate)] / sample rate
% ELSE
%   SET true start, duration = 0
% END
if(rate > 0)
    if (startWholeSec ~= gpsStart)
        rawFrac = gpsStart - startWholeSec;
        trueFrac = (floor(rawFrac*rate))/rate;
        trueStart = startWholeSec + trueFrac;
    else
        trueStart = startWholeSec;
    end
    trueDur = (floor(duration*rate))/rate;
else
    trueStart = 0;
    trueDur = 0;
end

return
