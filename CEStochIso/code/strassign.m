function d = strassign(in)
% STRASSIGN - convert argument to numeric if possible
%       
% input         output
%   5           5 (numeric)
%  '5'          5 (numeric)
%  'A'          A (string)
%  [1:5]        [1:5] (numeric vector)
%  '[1:5]'      [1:5] (numeric vector)
%  {'A1','K1','L2'} {'A1','K1','L2'} (cell array)
%  '[A1,K1,L2]' {'A1','K1','L2'} (cell array of strings
%  '[]'         empty array  
%  []           empty array 
%
% $Id: strassign.m,v 1.3 2007-02-27 16:16:35 kathorne Exp $

% IF input is a string
%   TRY to convert string to a number
%   IF result is numeric
%       SET output to converted result
%   ELSEIF string has vector format '[****]'
%       IF just a pair a brackets
%           SET output to MATLAB empty array []
%           EXIT
%       ENDIF
%       FIND location of separators (commas)
%       SET number of elements = # of commas + 1
%       CREATE cell array
%       LOOP over elements
%          FIND location of separators for element
%          GET sub-string between separators
%          SET cell element = sub-string
%       ENDLOOP
%   ELSE        
%       SET output to input
%   ENDIF
% ELSE
%   SET output to input
% ENDIF
if (ischar(in) && ~iscell(in))
    numVal = str2num(in);
    if(~isempty(numVal) && isnumeric(numVal))
        d = numVal;
    elseif ( (strcmp(in(1),'[')) && (strcmp(in(end),']')) )
        if(strcmp(in,'[]')==true)
            d = [];
            return
        end
        commaPos = strfind(in,',');
        endIn = length(in);
        numCell = numel(commaPos) + 1;
        d = cell(1,numCell);
        for icell = 1:numCell
            if (icell == 1)
                pStart = 1;
                if(numCell > 1)
                    pEnd = commaPos(1);
                else
                    pEnd = endIn;
                end
            elseif (icell == numCell)
                if(numCell > 1)
                    pStart = commaPos(end);
                else
                    pStart = 1;
                end
                pEnd = endIn;
            else
                pStart = commaPos(icell-1);
                pEnd = commaPos(icell);
            end
            cellStr = in((pStart+1):(pEnd-1));
            d{1,icell} = cellStr;
        end
    else
        d = cell(1,1);
        d{1,1} = in;
    end
else
    d = in;
end
return
