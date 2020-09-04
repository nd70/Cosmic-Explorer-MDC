function detector = buildbardetector(loc, axis);
%  BUILDBARDETECTOR -- Build Cartesian structure describing bar geometry
%
%  buildbardetector(loc, axis) builds a detector structure
%  describing the Cartesian location and response of an resonant bar
%  gravitational wave detector, given its location in geodetic
%  coordinates (in the structure loc) and the azimuth and altitude
%  angles of its long axis (in the structure axis).
%
%  The output is in the form of a structure with the fields
%      r: [3x1 double] %  position vector (in units of meters)
%                         in Earth-based Cartesian coordinates
%      d: [3x3 double] %  response tensor in Earth-based Cartesian coordinates
%
%  The input location is a structure loc with fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%             Meridian) in radians
%     height: elevation in meters above the WGS-84 reference ellipsoid 
%  Such a structure can be created from the geographic coordinates (in
%  degrees) using the function CREATELOCATION
%
%  The input axis direction is a structure axis with the fields
%        az: azimuth in radians East (clockwise) of North
%       alt: altitude (tilt) angle in radians above the local tangent plane
%  Such a structure can be created from the local angles (in degrees)
%  using the function CREATEORIENTATION
%  
%  The function calls GETCARTESIANPOSITION to determine the position
%  vector and GETBARRESPONSE to determine the response tensor and
%  then packs the results into a structure.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  See also GETCARTESIANPOSITION, GETBARRESPONSE, CREATELOCATION,
%  CREATEORIENTATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = getcartesianposition(loc);
d = getbarresponse(loc,axis);

detector = struct('r',r,'d',d);
