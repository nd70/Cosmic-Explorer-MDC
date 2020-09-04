function d = getiforesponse(loc,xarm,yarm);
%  GETIFORESPONSE -- Construct Cartesian response tensor from
%                    geographic coordinates and IFO arm orientations
%
%  getiforesponse(loc,xarm,yarm) calculates the Cartesian response
%  tensor associated with an interferometric gravitational wave detector
%  at a location loc with arms whose orientations are described by the
%  structures xarm and yarm.
%
%  The input location is a structure loc with fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%            Meridian) in radians
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
%  The function calls GETCARTESIANDIRECTION to convert the
%  orientation angles for each arm into a Cartesian unit vector
%  pointing along the arm, then constructs the response tensor as
%  one-half the difference between the outer products of these 
%  unit vectors with themselves.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
% 
%  See also GETCARTESIANDIRECTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = getcartesiandirection(xarm,loc);
v = getcartesiandirection(yarm,loc);

d = (u * transpose(u) - v * transpose(v) )./2;
