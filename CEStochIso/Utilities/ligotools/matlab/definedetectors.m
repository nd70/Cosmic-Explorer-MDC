%  DEFINEDETECTORS -- Define detector geometry structures for major detectors
%
%   This script defines structures containing geometry information for
%   the five IGEC resonant bar detectors and the five "first-generation"
%   interferometer sites:
%  
%   IFOs:
%     detector_LLO      % LIGO Livingston (Livingston, LA, USA)
%     detector_LHO      % LIGO Hanford (Hanford, WA, USA
%     detector_VIRGO    % VIRGO (Cascina, Italy)
%     detector_GEO600   % GEO-600 (Ruthe, Germany)
%     detector_TAMA     % TAMA-300 (Tokyo, Japan)
%  
%   Bars:
%     detector_ALLEGRO  % ALLEGRO (Baton Rouge, LA, USA)
%     detector_AURIGA   % AURIGA (Legnaro, Italy)
%     detector_EXPLORER % EXPLORER (Geneva, Switzerland)
%     detector_NAUTILUS % NAUTILUS (Frascati, Italy)
%     detector_NIOBE    % NIOBE (Perth, Australia)
%  
%   The fields of the IFO structures are
%       loc: [1x1 struct]  % Geographical location as generated by
%                            CREATELOCATION
%      xarm: [1x1 struct]  % X arm orientation as generated by
%                            CREATEORIENTATION
%      yarm: [1x1 struct]  % Y arm orientation as generated by
%                             CREATEORIENTATION
%       det: [1x1 struct]  % Detector geometry structure as generated
%                            by buildifodetector
%   See the referenced functions for more details
%  
%   The fields of the Bar structures are
%           loc: [1x1 struct]  % Geographical location as generated 
%                                by CREATELOCATION
%      IGECaxis: [1x1 struct]  % Orientation of long axis when roughly
%                                parallel to other IGEC bars, as 
%                                generated by CREATEORIENTATION
%       IGECdet: [1x1 struct]  % Detector geometry structure (in IGEC 
%                                orientation) as generated by buildbardetector
%   See the referenced functions for more details
%
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIGO Livingston (Livingston, Louisiana, USA)

detector_LLO.loc  = createlocation( 30+(33+46.4196/60)/60 , ...
				   - (90+(46+27.2654/60)/60) , ...
				   -6.574 );
detector_LLO.xarm = createorientation(180+72.2835, -3.121e-4*180/pi);
detector_LLO.yarm = createorientation(180-17.7165, -6.107e-4*180/pi);
detector_LLO.det  = buildifodetector(detector_LLO.loc, ...
                                      detector_LLO.xarm, ...
                                      detector_LLO.yarm);

% LIGO Hanford (Hanford, Washington, USA)

detector_LHO.loc  = createlocation( 46+(27+18.528/60)/60 , ...
				   - (119+(24+27.5657/60)/60) , ...
				   142.554);
detector_LHO.xarm = createorientation(-35.9994, -6.195e-4*180/pi);
detector_LHO.yarm = createorientation(180+54.0006, 1.25e-5*180/pi);
detector_LHO.det  = buildifodetector(detector_LHO.loc, ...
                                      detector_LHO.xarm, ...
                                      detector_LHO.yarm);

% GEO-600 (Hannover, Germany)

detector_GEO600.loc  = createlocation(52 + (14 + 42.528/60)/60 , ...
				      9 + (48 + 25.894/60)/60 , ...
				      114.425);
detector_GEO600.xarm = createorientation(90-21.6117);
detector_GEO600.yarm = createorientation(90-115.9431);
detector_GEO600.det  = buildifodetector(detector_GEO600.loc, ...
                                      detector_GEO600.xarm, ...
                                      detector_GEO600.yarm);

% TAMA-300 (Tokyo, Japan)

detector_TAMA300.loc  = createlocation(35 + (40 + 35.6/60)/60 , ...
				       139 + (32 + 9.8/60)/60 , ...
				       90);
detector_TAMA300.xarm = createorientation(90-180);
detector_TAMA300.yarm = createorientation(90-270);
detector_TAMA300.det  = buildifodetector(detector_TAMA300.loc, ...
                                      detector_TAMA300.xarm, ...
                                      detector_TAMA300.yarm);

% VIRGO (Cascina/Pisa, Italy)

detector_VIRGO.loc  = createlocation(43 + (37 + 53.0921/60)/60 , ...
				     10 + (30 + 16.1878/60)/60 , ...
				     51.884);
detector_VIRGO.xarm = createorientation(90-70.5674);
detector_VIRGO.yarm = createorientation(90-160.5674);
detector_VIRGO.det  = buildifodetector(detector_VIRGO.loc, ...
                                      detector_VIRGO.xarm, ...
                                      detector_VIRGO.yarm);

% ALLEGRO (Baton Rouge, Louisiana, USA)
detector_ALLEGRO.loc = createlocation( 30+(24+45.110/60)/60, ...
				      - (91+(10+43.766/60)/60) );
detector_ALLEGRO.IGECaxis = createorientation(-40);
detector_ALLEGRO.IGECdet  = buildbardetector(detector_ALLEGRO.loc, ...
					      detector_ALLEGRO.IGECaxis);

% AURIGA (Legnaro/Padova, Italy)
detector_AURIGA.loc = createlocation( 45+(21+12/60)/60, ...
                                     11+(56+54/60)/60, ...
                                     0 );
detector_AURIGA.IGECaxis = createorientation(44);
detector_AURIGA.IGECdet  = buildbardetector(detector_AURIGA.loc, ...
					     detector_AURIGA.IGECaxis);

% EXPLORER (Geneva, Switzerland)
detector_EXPLORER.loc = createlocation( 46+27/60, ...
				       6+12/60 );
detector_EXPLORER.IGECaxis = createorientation(39);
detector_EXPLORER.IGECdet  = buildbardetector(detector_EXPLORER.loc, ...
					       detector_EXPLORER.IGECaxis);

% NAUTILUS (Frascati/Rome, Italy)
detector_NAUTILUS.loc = createlocation( 41+(49+26/60)/60, ...
				       12+(40+21/60)/60 );
detector_NAUTILUS.IGECaxis = createorientation(44);
detector_NAUTILUS.IGECdet  = buildbardetector(detector_NAUTILUS.loc, ...
					    detector_NAUTILUS.IGECaxis);

% NIOBE (Perth, Australia)
detector_NIOBE.loc = createlocation(-(31+56/60), ...
				       115+49/60 );
detector_NIOBE.IGECaxis = createorientation(0);
detector_NIOBE.IGECdet  = buildbardetector(detector_NIOBE.loc, ...
					       detector_NIOBE.IGECaxis);
