function loc = createlocation(latDeg, lonDeg, h);
%  CREATELOCATION -- Create geodetic location structure
%                    (converting angles from degrees to radians)
%
%  createlocation(latDeg, lonDeg, h) creates a structure with the fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%             Meridian) in radians
%     height: elevation in meters above the WGS-84 reference ellipsoid 
%  from the latitude and longitude (in degrees) and the height in
%  meters above the reference ellipsoid.  If the third argument is
%  omitted, the height is set to zero.
%
%  This function simply converts the input angles from degrees into
%  radians, sets h to zero if it's not included in the argument list,
%  and then packs the latitude, longitude and height into a structure.
%
%  Routine written by John T. Whelan.
%  Contact jtwhelan@loyno.edu
%  $Id: createlocation.m,v 1.2 2004/01/31 18:55:17 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2)
  h = 0;
end
loc = struct('lat',latDeg*pi/180,'lon',lonDeg*pi/180,'height',h);
