function site = getsitefromletter(letter)
% GETSITEFROMLETTER -- get detector site string corresponding to
%                      first letter in channel prefix
%
% sitefromletter(letter) returns the detector site string corresponding to
% a given first letter of a channel prefix (e.g., 'H' in 'H1:LSC-AS_Q'
% corresponds to 'LHO').  The correspondence is as specified in the
% frame format specification,
% http://www.ligo.caltech.edu/docs/T/T970130-F.pdf
% 
% L = LLO
% H = LHO
% V = VIRGO
% G = GEO600
% T = TAMA
% A = ALLEGRO
% O = AURIGA
% E = EXPLORER
% N = NAUTILUS
% B = NIOBE
% 
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%  $Id: getsitefromletter.m,v 1.3 2006-04-14 21:05:42 whelan Exp $
%
%  See also GETDETECTOR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch letter

  case 'L'
    site = 'LLO';

  case 'H';
    site = 'LHO';

  case 'V';
    site = 'VIRGO';

  case 'G';
    site = 'GEO600';

  case 'T';
    site = 'TAMA';

  case 'A';
    site = 'ALLEGRO';

  case 'O';
    site = 'AURIGA';

  case 'E';
    site = 'EXPLORER';

  case 'N';
    site = 'NAUTILUS';

  case 'B';
    site = 'NIOBE';

  otherwise
    error(['invalid detector site letter ' letter]);

end;

return
