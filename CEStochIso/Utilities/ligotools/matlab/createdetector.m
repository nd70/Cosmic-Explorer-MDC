function detector = createdetector(r, d);
%  CREATEDETECTOR -- Create Cartesian detector structure from position
%                    and response
%
%  createdetector(r, d) builds a detector structure out of the
%  Cartesian position vector r and response tensor d
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

detector = struct('r',r,'d',d);
