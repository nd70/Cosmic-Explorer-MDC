function data = ilwdread(filename, verbose)
%ILWDREAD Load an ILWD (Internal LightWeight Data) formatted XML file.
%   ILWDREAD(FILENAME) returns a MATLAB struct representation of the ILWD file
%   given by FILENAME.  Each element of the ILWD is a MATLAB struct with field
%   "tag" storing the element type, and additional fields storing the values of
%   the attributes of the element.  If the element is an array (such as real_8)
%   a field named after the element ("real_8") stores the array data in MATLAB
%   format.  If the element is a container (such as ilwd) the member elements
%   are stored as fields m1, m2, ...
%
%   ILWDREAD(FILENAME, VERBOSE) turns on verbose output if VERBOSE is true
%   (nonzero).  This consists of a formatted representation of the structure
%   but not the content of the ILWD (i.e. the tags and fields but not the
%   data).
%
%   The structure of ILWD files is often complex and this is reflected in the
%   complicated structs produced by ILWDREAD.  The struct can be navigated by
%   typing the element names on the command line to see the available fields.
%
%   >>>> ilwd = ilwdread('results.ilwd')
%
%   ilwd = 
%   
%           tag: 'ilwd'
%       comment: 'ilwdread example'
%          name: 'ligo:ldas:file'
%            m1: [1x1   struct]
%   
%   >>>> ilwd.m1
%   
%   ans = 
%   
%          tag: 'int_4s'
%         dims: 2
%         name: 'one two'
%       int_4s: [2x1 double]
%
%   >>>> ilwd.m1.int_4s
%
%   ans =
%
%        1
%        2
%   
%   It may be convenient to copy the data to a concisely named variable.
%
%   >>>> a = ilwd.m1.int_4s;
%
%   The member structs fields of a struct must be name-mangled (m1, m2, ...)
%   rather than subscripted (m(1), m(2), ...) to allow them to represent
%   different structures.  The EVAL command can be used to simulate
%   subscripting, replacing "data.m(i)" with "eval(['data.m' i])".
%
%   The array fields are MATLAB vectors containing numerical data.  The
%   values are preserved but the format is not (for example, int types are
%   converted to double).  Complex types are also unpacked to construct 
%   complex vectors.
%
%   ILWDREAD accepts only uncompressed ASCII ILWD files.  Any non-standard
%   attributes should be correctly parsed to a string (but will issue a
%   warning).  Any non-standard tags will cause an error.  Bugs and comments
%   to Antony.Searle@anu.edu.au
%
%   See also

%   Author : Antony C. Searle
%   Version: 0.2

fid = fopen(filename, 'rt');            % open the file
if fid == -1                            % check for errors
   error('open failed');
end
header = fscanf(fid, '%c', 8);          % read in the file header
if nargin < 2                           % check for a verbose argument
   verbose = 0;                         % no argument so set to default (quiet)
end
if header ~= '<?ilwd?>'                 % confirm we have an ilwd file
   error(['unhandled format ' header]);
end
if verbose
   disp(header);
end
data = parsetag(fid, '', verbose);      % recursively parse the ilwd element
fclose(fid);                            % close the file

% parsetag helper function

function data = parsetag(fid, indent, verbose);
%PARSETAG Recursively parse ILWD tags for ILWDREAD
%   See also ILWDREAD

%   Author : Antony C. Searle
%   Version: 0.2

c = fscanf(fid, '%c', 1);               % find the tag start
while c ~= '<'
   c = fscanf(fid, '%c', 1);
end

[data.tag, terminator] = parseidentifier(fid);
                                        % read the tag

local_ndim = 1;                         % set default attribute values
local_size = 1;
local_dims = 1;

if verbose
   disp([indent '<' data.tag]);
end

flag = 1;                                   
while flag                              % parse attributes
   if terminator ~= '>'                 % if not at the tag end, parse an
                                        %    attribute
      [string, terminator] = parseidentifier(fid);
   end
   switch terminator                    % switch on identifier terminating
                                        %    character
   case '='                             % assigned attribute
      fscanf(fid, '%c', 1);             % ignore opening quote
      switch string                     % switch on attribute name
      case 'size'
         data.size = fscanf(fid, '%d', 1);
                                        % read decimal integer
         local_size = data.size;
         fscanf(fid, '%c', 1);          % ignore the closing quote
         if verbose
            disp([indent '    size=' num2str(data.size)]);
         end
      case 'ndim'
         data.ndim = fscanf(fid, '%d', 1); 
                                        % read decimal integer
         local_ndim = data.ndim;
         fscanf(fid, '%c', 1);          % ignore the closing quote
         if verbose
            disp([indent '    ndim=' num2str(data.ndim)]);
         end
      case 'dims'
         data.dims = 1;                 % touch the field
         for i = 1:local_ndim           % loop over dimensions
            data.dims(i) = fscanf(fid, '%d', 1);
                                        % read in dimensions
            fscanf(fid, '%c', 1);       % ignore the closing quote
         end
         local_dims = data.dims;
         if verbose
            disp([indent '    dims=' num2str(data.dims)]);
         end
      case 'comment'
         [data.comment, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    comment=' data.comment]);
         end
      case 'name'
         [data.name, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    name=' data.name]);
         end
      case 'units'
         [data.units, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    units=' data.units]);
         end
      case 'mdorder'
         [data.mdorder, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    mdorder=' data.mdorder]);
         end
      case 'nullmask'
         [data.nullmask, terminator] = parseattribute(fid);
                                        %read string
         if verbose
            disp([indent '    nullmask=' data.nullmask]);
         end
      case 'format'
         [data.format, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    format=' data.format]);
         end
         if data.format ~= 'ascii'      % check format is ASCII
            error(['unhandled format: ' data.format]);
         end
      case 'compression'
         [data.compression, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    format=' data.format]);
         end
         if (data.compression ~= 'none') | (data.compression ~= 'gzip0');
                                        % check data is uncompressed
            error(['unhandled compression: ' data.compression']);
         end
      case 'byteorder'
         error(['unhandled binary attribute: byteorder']);
                                        % implies unsupported binary format
      case 'bytes'
         error(['unhandled binary attribute: bytes']);
                                        % implies unsupported binary format
      case 'parser'
         [data.parser, terminator] = parseattribute(fid);
                                        % read string
         if verbose
            disp([indent '    parser=' data.parser]);
         end
      case 'metadata'
          [data.metadata, terminator] = parseattribute(fid);
                                        % read string
          if verbose
              disp([indent '    metadata=' data.metadata]);
          end
      case 'dataValueUnit'
          [data.dataValueUnit, terminator] = parseattribute(fid);
                                        % read string
          if verbose
              disp([indent '    dataValueUnit=' data.dataValueUnit]);
          end
      case 'dx'
         data.dx = fscanf(fid, '%f', 1);
                                        % read floating-point number
         fscanf(fid, '%c', 1);          % ignore the closing quote
         if verbose
            disp([indent '    dx=' num2str(data.dx)]);
         end         
     case 'startx'
         data.startx = fscanf(fid, '%f', 1);
                                        % read floating-point number
         fscanf(fid, '%c', 1);          % ignore the closing quote
         if verbose
            disp([indent '    startx=' num2str(data.startx)]);
         end
      otherwise
         disp(['warning: nonstandard attribute ' string]);
         eval(['[data.' string ', terminator] = parseattribute(fid);']);
                                        % read string into nonstandard attribute
      end
   case '>'                             % end of the tag
      if verbose
         disp([indent '    >']);
      end
      flag = 0;                         % stop while loop
   otherwise
      error(['unhandled character ' terminator]);
                                        % unexpected character
   end
end

switch data.tag                         % switch on tag type
case 'ilwd'                             % container type
   indent2 = ['    ' indent];           % set indentation level
   for i = 1:local_size                 % read in "size" elements
      eval(['data.m' num2str(i) ' = parsetag(fid, indent2, verbose);']);
                                        % construct the field name and parse
                                        %    into it
   end
case 'char'
   data.char = fscanf(fid, '%c', data.dims);
                                        % read in "dims" characters
case 'char_u'
   buffer = fscanf(fid, '%c%d');        % read in as many character-separated
                                        % decimal integers as possible
   data.char_u = buffer(2:2:length(buffer));
                                        % strip out the integers                                      
case 'char_s'
   buffer = fscanf(fid, '%c%d');        % read in as many character-separated
                                        % decimal integers as possible
   data.char_s = buffer(2:2:length(buffer));
                                        % strip out the integers                                      
case 'int_2u'
   data.int_2u = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'int_4u'
   data.int_4u = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'int_8u'
   data.int_8u = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'int_2s'
   data.int_2s = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'int_4s'
   data.int_4s = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'int_8s'
   data.int_8s = fscanf(fid, '%d');     % read in as many decimal integers as
                                        % possible
case 'real_4'
   data.real_4 = fscanf(fid, '%f');     % read in as many decimal integers as
                                        % possible
case 'real_8'
   data.real_8 = fscanf(fid, '%f');     % read in as many decimal integers as
                                        % possible
case 'complex_8'
   buffer = fscanf(fid, '%f');          % read in as many floating point 
                                        %    decimals as possible
   data.complex_8 = buffer(1:2:size(buffer)) + buffer(2:2:size(buffer)) * ...
      sqrt(-1);                         % unpack to a complex array
case 'complex_16'
   buffer = fscanf(fid, '%f');          % read in as many floating point 
                                        %    decimals as possible
   data.complex_16 = buffer(1:2:size(buffer)) + buffer(2:2:size(buffer)) * ...
      sqrt(-1);                         % unpack to a complex array
case 'lstring'
   data.lstring = fscanf(fid, '%c', local_size);
                                        % read "size" characters into a string
case 'external'
   data.external = fscanf(fid, '%c', local_size);
                                        % read "size" characters into a string
otherwise
   error(['unrecognised tag ' data.tag]);
                                        % tag is not supported
end

[string, terminator] = parseidentifier(fid);      
                                        % ignore the end tag
if verbose
   disp([indent string terminator]);
end

% parseattribute helper function

function [string, terminator] = parseattribute(fid);
%PARSEATTRIBUTE Parse attributes for ILWDREAD
%   See also ILWDREAD

%   Author : Antony C. Searle
%   Version: 0.2

string = '';                            % empty output
flag = 1;
while flag;                             % read in characters
   terminator = fscanf(fid, '%c', 1);   % get a character
   switch terminator                    % switch on the character
   case double(39)                      % single-quote terminates
      flag = 0;
   case double(34)                      % double-quote terminates
      flag = 0;
   otherwise                            % any other character
      string = [string terminator];     %    is part of output
   end
end

% parseidentifer

function [string, terminator] = parseidentifier(fid);
%PARSEIDENTIFIER Parse identifiers for ILWDREAD
%   See also ILWDREAD

%   Author : Antony C. Searle
%   Version: 0.2

string = '';                            % empty output
flag = 1;
white = 1;                              % leading whitespace?
while flag;                             % read in characters
   terminator = fscanf(fid, '%c', 1);   % get a character
   switch terminator                    % switch on the character
   case '='                             % equals terminates
      flag = 0;
   case '>'                             % right-angle terminates
      flag = 0;
   case double(39)                      % single-quote terminates
      flag = 0;
   case double(34)                      % double-quote terminates
      flag = 0;
   case double(10)                      % newline terminates unless
                                        %    still leading whitespace
      flag = white;
   case ' '                             % space terminates unless
                                        %    still leading whitespace
      flag = white;
   otherwise                            % any other character
      white = 0;                        %    terminates whitespace
      string = [string terminator];     %    and is part of output
   end
end







