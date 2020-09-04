function detector = buildifodetector(loc, xarm, yarm);
%  BUILDIFODETECTOR -- Build Cartesian structure describing IFO geometry
%
%  buildifodetector(loc, xarm, yarm) builds a detector structure
%  describing the Cartesian location and response of an interferometric
%  gravitational wave detector, given its location in geodetic
%  coordinates (in the structure loc) and the  azimuth and altitude
%  angles of its arms (in the structures xarm and yarm).
%
%  The output is in the form of a structure with the fields
%      r: [3x1 double] % position vector (in units of meters)
%                        in Earth-based Cartesian coordinates
%      d: [3x3 double] % response tensor in Earth-based Cartesian coordinates
%
%  The input location is a structure loc with fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%             Meridian) in radians
%     height: elevation in meters above the WGS-84 reference ellipsoid 
%  Such a structure can be created from the geographic coordinates (in
%  degrees) using the function CREATELOCATION
%
%  The input arm directions are structures xarm and yarm with the fields
%        az: azimuth in radians East (clockwise) of North
%       alt: altitude (tilt) angle in radians above the local tangent plane
%  Such a structure can be created from the local angles (in degrees)
%  using the function CREATEORIENTATION
% 
%  The function calls GETCARTESIANPOSITION to determine the position
%  vector and GETIFORESPONSE to determine the response tensor and
%  then packs the results into a structure.
%
%  Routine written by John T. Whelan.
%  Contact jtwhelan@loyno.edu
%  $Id: buildifodetector.m,v 1.3 2004/01/31 18:55:17 whelan Exp $
% 
%  See also GETCARTESIANPOSITION, GETIFORESPONSE, CREATELOCATION,
%  CREATEORIENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = getcartesianposition(loc);
d = getiforesponse(loc,xarm,yarm);

detector = struct('r',r,'d',d);
