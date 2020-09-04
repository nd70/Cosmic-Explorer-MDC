function initCombine(params)

% initializes the variable for live combining of results of isotropic analysis
% Previously done in the old stochastic.m
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: initCombine.m,v 1.2 2007-07-20 22:44:16 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global COMB;

if params.doCombine
  %initialize various variables
  if params.doOverlap
    COMB.add_ccstats_o = 0;
    COMB.add_ccspectra_o = 0;
    COMB.add_sensint_o = 0;
    COMB.add_ccstats_e = 0;
    COMB.add_ccspectra_e = 0;
    COMB.add_sensint_e = 0;
    COMB.dataI = 0;
    COMB.dataJ = 0;
    COMB.addctr = 1;
    COMB.ccVars = zeros(1,params.numIntervalsTotal);
    COMB.combined_ccstats = 0;
    COMB.combined_ccspectra = 0;
    COMB.combined_sensint = 0;
    COMB.ccVars_ovl = zeros(1,params.numIntervalsTotal);
    COMB.ovladdctr = 1;
  else
    COMB.add_ccstats = 0;
    COMB.add_ccspectra = 0;
    COMB.add_sensint = 0;
    COMB.addctr = 1;
    COMB.ccVars = zeros(1,params.numIntervalsTotal);
  end

  if params.doBadGPSTimes
    COMB.badGPSTimesctr = 1;
    COMB.badGPSTimes = zeros(1,params.numIntervalsTotal);
    if 1
      if isempty(params.badGPSTimesFile)
        COMB.badtimesstart = 9999999999;
        COMB.badtimesend = 0;
      else
        [COMB.badtimesstart,COMB.badtimesend] = textread(params.badGPSTimesFile,'%f%f\n',-1,'commentstyle','matlab');
      end
    else
      if isempty(params.badGPSTimesFile)
        COMB.badtimesstart = 9999999999;
        COMB.badtimesend = 0;
      else
        COMB.badtimesstart = textread(params.badGPSTimesFile,'%f\n',-1,'commentstyle','matlab');
        COMB.badtimesend = COMB.badtimesstart+60;
      end
    end
  end
end
