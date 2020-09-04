function d = getbarresponse(loc,axis);
%  GETBARRESPONSE -- Construct Cartesian response tensor from
%                    geographic coordinates & bar orientation
%
%  getbarresponse(loc,axis) calculates the Cartesian response tensor
%  associated with a resonant bar gravitational wave detector at a
%  location loc with arms whose long axis has an orientation described
%  by the structure's axis.
%
%  The input location is a structure loc with fields
%        lat: geodetic latitude (measured North from the Equator) in radians
%        lon: geodetic longitude (measured East from the Prime
%            Meridian) in radians
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
%  The function calls GETCARTESIANDIRECTION to convert the
%  orientation angles for the axis into a Cartesian unit vector
%  pointing along the axis, then constructs the response tensor as
%  the outer product of this unit vector with itself.
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
% 
%  See also GETCARTESIANDIRECTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = getcartesiandirection(axis,loc);

d = u * transpose(u);
