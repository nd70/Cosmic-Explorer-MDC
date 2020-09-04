function detector = getdetector(site, azDeg)
% GETDETECTOR -- get detector geometry structure for major detector
%
% getdetector(site, azDeg) returns the detector structure for the detector
% referred to by the string site.  Recognized values for site are
% 
% LLO, LHO, VIRGO, GEO600, TAMA,             % (IFOs)
% ALLEGRO, AURIGA, EXPLORER, NAUTILUS, NIOBE % (bar detectors)
%
% The optional argument azDeg (only allowed for bar detectors)
% specifies the orientation of the bar in degrees East (clockwise)
% from local North.  If it is omitted, the orientation roughly 
% parallel to the other IGEC detectors is used.
%
%  The output is in the form of a structure with the fields
%      r: [3x1 double] %  position vector (in units of meters)
%                         in Earth-based Cartesian coordinates
%      d: [3x3 double] %  response tensor in Earth-based Cartesian coordinates
%
% The function calls BUILDIFODETECTOR or BUILDBARDETECTOR as appropriate
% to construct the detector structure using hard-coded geographical data.
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  See also BUILDIFODETECTOR, BUILDBARDETECTOR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch site

  case 'LLO'
    % LIGO Livinston (Livingston, Louisiana, USA)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation( 30+(33+46.4196/60)/60 , ...
			  - (90+(46+27.2654/60)/60) , ...
   			  -6.574 );
    xarm = createorientation(180+72.2835, -3.121e-4*180/pi);
    yarm = createorientation(180-17.7165, -6.107e-4*180/pi);
    detector = buildifodetector(loc, xarm, yarm);

  case 'LHO'
    % LIGO Hanford (Hanford, Washington, USA)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation( 46+(27+18.528/60)/60 , ...
			  - (119+(24+27.5657/60)/60) , ...
			  142.554);
    xarm = createorientation(-35.9994, -6.195e-4*180/pi);
    yarm = createorientation(180+54.0006, 1.25e-5*180/pi);
    detector = buildifodetector(loc, xarm, yarm);

  case 'VIRGO'
    % VIRGO (Cascina/Pisa, Italy)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation(43 + (37 + 53.0921/60)/60 , ...
			  10 + (30 + 16.1878/60)/60 , ...
			  51.884);
    xarm = createorientation(90-70.5674);
    yarm = createorientation(90-160.5674);
    detector = buildifodetector(loc, xarm, yarm);

  case 'GEO600'
    % GEO-600 (Hannover, Germany)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation(52 + (14 + 42.528/60)/60 , ...
			  9 + (48 + 25.894/60)/60 , ...
			  114.425);
    xarm = createorientation(90-21.6117);
    yarm = createorientation(90-115.9431);
    detector = buildifodetector(loc, xarm, yarm);

  case 'TAMA300'
    % TAMA-300 (Tokyo, Japan)

    if (nargin ~= 1)
      error('Cannot specify orientation of IFO\n');
    end

    loc  = createlocation(35 + (40 + 35.6/60)/60 , ...
			  139 + (32 + 9.8/60)/60 , ...
			  90);
    xarm = createorientation(90-180);
    yarm = createorientation(90-270);
    detector = buildifodetector(loc, xarm, yarm);

  case 'ALLEGRO'
    % ALLEGRO (Baton Rouge, Louisiana, USA)

    if (nargin == 1)
      azDeg = -40;
    end

    loc = createlocation(30+(24+45.110/60)/60, ...
			 - (91+(10+43.766/60)/60) );
    detector = buildbardetector( loc, createorientation(azDeg) );

  case 'AURIGA'
    % AURIGA (Legnaro/Padova, Italy)

    if (nargin == 1)
      azDeg = 44;
    end

    loc = createlocation(45+(21+12/60)/60, ...
                         11+(56+54/60)/60, ...
                         0 );
    detector = buildbardetector( loc, createorientation(azDeg) );

  case 'EXPLORER'
    % EXPLORER (Geneva, Switzerland)

    if (nargin == 1)
      azDeg = 39;
    end

    loc = createlocation(46+27/60, ...
			 6+12/60 );
    detector = buildbardetector( loc, createorientation(azDeg) );

  case 'NAUTILUS'
    % NAUTILUS (Frascati/Rome, Italy)

    if (nargin == 1)
      azDeg = 44;
    end

    loc = createlocation(41+(49+26/60)/60, ...
			 12+(40+21/60)/60 );
    detector = buildbardetector( loc, createorientation(azDeg) );

  case 'NIOBE'

    if (nargin == 1)
      azDeg = 0;
    end

    % NIOBE (Perth, Australia)
    loc = createlocation(-(31+56/60), ...
			 115+49/60 );
    detector = buildbardetector( loc, createorientation(azDeg) );

  otherwise

    try
      detector = getdetector(getsitefromletter(site),azDeg);
    catch
      try
	detector = getdetector(getsitefromletter(site));
      catch
	error(['invalid detector site ' site]);
      end;
    end;

end;

return;
