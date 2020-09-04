function chStruct = chanstruct(chanName,startTime)
%CHANSTRUCT - Channel structure contructor
%  chStruct = chanstruct(chanName[,startTime]) creates an IGWD Channel
%                   structure using the channel name 'chanName'
%      The syntax of valid channel names (i.e. H2:LSC-AS_Q)
%     Xn:Name X is IFO location (H Hanford, L Livingston, V Virgo)
%             n is detector number (0 for environment monitoring)
%             Name is detector channel name, usually a location and
%                   signal type
%  startTime - [OPTIONAL] Sets a base GPS start time (gpsStart) so that
%    all time-series ranges are offset from this base (instead of base = 0)
%
%  chStruct has the following structure
%       name - channel name (from chanName input)
%       statusCode - 0 if name syntax OK, 1 otherwise
%       site - data location (default LDR) **UNUSED AT PRESENT**
%       instrument - IFO site (default first letter of channel name)
%       type -  frame data type (default RDS_R_L1)
%       gpsStart - base GPS start time (default = 0)
%       rate - channel sample rate (default = 0)
% $Id: chanstruct.m,v 1.1 2004-08-11 18:43:34 kathorne Exp $

% IF invalid # of inputs, outputs
%   flag as error
%   exit
% ENDIF
inErr = nargchk(1,2,nargin);
outErr = nargoutchk(1,1,nargout);
if(~isempty(inErr) || ~isempty(outErr))
  msgid = 'chanstruct:inOutArg';
  errmsg = sprintf('%s:\n\t%s',msgid,...
      'chStruct=chanstruct(chanName,[startTime])');
  error(msgid,errmsg);
end

% create Channel structure
% SET name from first argument
% IF name has correct syntax
%   SET status code = 0
% ELSE
%   SET status code = 1
% ENDIF
% SET default site
% SET default instrument, type
% IF there is a second argument
%   IF the startTime is valid (> 0)
%       SET start time to arguement value
%   ELSE
%       SET default start time
%   ENDIF
% ELSE
%   SET default start time
% ENDIF
% SET rate = 0
chStruct = struct('name',' ','statusCode',1,'site',' ','instrument',' ',...
        'type',' ','gpsStart',0,'rate',0);
chStruct.name = chanName;
if(isletter(chanName(1:1)) & ~isletter(chanName(2:2)) & chanName(3:3) == ':')
    chStruct.statusCode = 0;
else
    chStruct.statusCode = 1;
end
chStruct.site = 'LDR';
chStruct.instrument = chanName(1);
chStruct.type = 'RDS_R_L1';
if (nargin == 2)
    if (startTime > 0)
        chStruct.gpsStart = startTime;
    else
        chStruct.gpsStart = 0;
    end
else
    chStruct.gpsStart = 0;
end
chStruct.rate = 0;
return
