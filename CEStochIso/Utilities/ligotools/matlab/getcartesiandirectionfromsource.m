function n = getcartesiandirectionfromsource(decDeg,minusHAHr);
%  GETCARTESIANDIRECTIONFROMSOURCE -- Construct Cartesian unit vector 
%                                     for wave propagation from
%                                     declination and minus hour angle
%                                  
%  
%  getcartesiandirectionfromsource(decDeg,minusHAHr) calculates the unit
%  vector n in Earth-based Cartesian coordinates corresponding to a
%  given declination and minus hour angle (which is equal to the
%  source's right ascension minus the Greenwich mean sidereal time
%  (GMST) at which it is observed).
%
%  The output is a unit vector u in an Earth-fixed Cartesian
%  coordinate system whose first, second, and third axes pierce the
%  surface of the Earth at 
%         1. The intersection of the Equator and the Prime Meridian
%         2. The intersection of the Equator and the meridian at 90
%            degrees East longitude
%         3. The North Pole
%
%  The inputs are
%     decDeg: declination in degrees north of the celestial equator
%  minusHAHr: minus the hour angle (i.e., right ascension minus GMST)
%             in hours EAST of Greenwich
%  
%  The function first converts the angles into radians and then
%  calculates the components of the unit vector pointed from the
%  source to the center of the Earth-based coordinate system
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  See also GETCARTESIANDIRECTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cosdelta = cos(decDeg*pi/180);
  
n = [-cosdelta.*cos(minusHAHr*pi/12);...
     -cosdelta.*sin(minusHAHr*pi/12);...
     -sin(decDeg*pi/180)];
return;
  