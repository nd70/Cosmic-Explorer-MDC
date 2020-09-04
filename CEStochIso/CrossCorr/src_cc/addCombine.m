function addCombine(params, I, segmentStartTime, result)

% adds the current result to the stored stuff
% Previously done in the old stochastic.m
% 
% Inputs:
%   params - structure containing the parameters used in the pipeline
%   I - the index of the current interval
%   segmentStartTime - the start time of the segment
%   result - structure containing results of pipeline (ccStat etc)
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: addCombine.m,v 1.2 2007-08-03 21:43:53 vmandic Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The COMB variable is global to initCombine, addCombine and finishCombineto
% facilitate data combination. It would be nice to find a more elegant solution.
global COMB;

%SHIFT: have to worry about different trials
    if params.doCombine
      if params.doOverlap

        if params.doBadGPSTimes %determine if the current interval is bad
          if 1
           c1 = params.intervalStartTime - params.bufferSecsMax < COMB.badtimesend;
           c2 = params.intervalStartTime + 3*params.segmentDuration + params.bufferSecsMax > COMB.badtimesstart;
          else
           c1 = params.intervalStartTime+60 == COMB.badtimesstart;
           c2 = 1;
          end

          if sum(c1&c2) > 0 | sqrt(result.ccVar/result.naiVar)>params.maxDSigRatio | sqrt(result.ccVar/result.naiVar)<params.minDSigRatio
            COMB.badtimesflag = 1;
            COMB.badGPSTimes(COMB.badGPSTimesctr) = segmentStartTime;
            COMB.badGPSTimesctr = COMB.badGPSTimesctr + 1;
          else
            COMB.badtimesflag = 0;
          end
        else
          COMB.badtimesflag = 0;
        end

	if ~COMB.badtimesflag %the interval is not bad, so add it
	  % crowder: expanded conditional statement to accomodate doDirectional being true or false
          if mod(COMB.addctr,2) == 1 & params.doDirectional == 1 %odd interval in this overlapping section and doDirectional=true
            COMB.add_ccstats_o = COMB.add_ccstats_o + result.ccStat.data / result.ccVar / params.segmentDuration;
            COMB.add_ccspectra_o = COMB.add_ccspectra_o + result.ccSpec.data / result.ccVar / params.segmentDuration;
            COMB.add_sensint_o = COMB.add_sensint_o + result.sensInt.data * params.segmentDuration^2;
	  elseif mod(COMB.addctr,2) == 1 & params.doDirectional == 0 %odd interval in this overlapping section and doDirectional=false
            COMB.add_ccstats_o = COMB.add_ccstats_o + result.ccStat / result.ccVar / params.segmentDuration;
            COMB.add_ccspectra_o = COMB.add_ccspectra_o + result.ccSpec.data / result.ccVar /params.segmentDuration;
            COMB.add_sensint_o = COMB.add_sensint_o + result.sensInt.data * params.segmentDuration^2;
          elseif mod(COMB.addctr,2) ~= 1 & params.doDirectional == 1 %even interval in this job and doDirectional=true
            COMB.add_ccstats_e = COMB.add_ccstats_e + result.ccStat.data / result.ccVar / params.segmentDuration;;
            COMB.add_ccspectra_e = COMB.add_ccspectra_e + result.ccSpec.data / result.ccVar /params.segmentDuration;
            COMB.add_sensint_e = COMB.add_sensint_e + result.sensInt.data * params.segmentDuration^2;
          else %even interval in this job and doDirectional=false
            COMB.add_ccstats_e = COMB.add_ccstats_e + result.ccStat / result.ccVar / params.segmentDuration;;
            COMB.add_ccspectra_e = COMB.add_ccspectra_e + result.ccSpec.data / result.ccVar / params.segmentDuration;
            COMB.add_sensint_e = COMB.add_sensint_e + result.sensInt.data * params.segmentDuration^2;
          end

          COMB.tmpsensintadd = result.sensInt.data * params.segmentDuration^2;
          COMB.dataI = COMB.dataI + COMB.tmpsensintadd;
          if COMB.addctr > 1
	    COMB.dataJ = COMB.dataJ + COMB.tmpsensintadd;
          end
          COMB.ccVars(COMB.addctr) = result.ccVar;
          COMB.addctr = COMB.addctr + 1;
        end

        if COMB.badtimesflag | I == params.numIntervalsTotal
	    %if the intervals is bad, or if we reached the last interval
	    %of the job, finish combining the last overlapping part of
	    %this job, and add it to the combined variables.

          COMB.ccVars(COMB.addctr:end) = [];
          if COMB.addctr>2
            if mod(COMB.addctr-1,2)==1
              j_o = [1:2:COMB.addctr-1]';
              j_e = [2:2:COMB.addctr-1-1]';
            else
              j_o = [1:2:COMB.addctr-1-1]';
              j_e = [2:2:COMB.addctr-1]';
            end

            errorBar_o = sqrt( 1/sum(1./COMB.ccVars(j_o)) ) / params.segmentDuration;
            errorBar_e = sqrt( 1/sum(1./COMB.ccVars(j_e)) ) / params.segmentDuration;

            C_oo = errorBar_o^2;
            C_ee = errorBar_e^2;
  
            [w2bar,w4bar,woverlap4bar]=windowFactors(params.fft1.dataWindow,params.fft2.dataWindow);
            sigma2I  = COMB.ccVars(1:end-1) / (params.segmentDuration^2);
            sigma2J  = COMB.ccVars(2:end)   / (params.segmentDuration^2);
            sigma2IJ = 0.5*(woverlap4bar/w4bar)*0.5*(sigma2I+sigma2J);
  
            C_oe = C_oo*C_ee*sum(sigma2IJ./(sigma2I.*sigma2J));
            detC = C_oo*C_ee - C_oe*C_oe;

            COMB.dataI = COMB.dataI - COMB.tmpsensintadd;
            COMB.dataIJ = 0.5*(woverlap4bar/w4bar)*0.5*(COMB.dataI + COMB.dataJ); 
            lambda_o = (C_ee - C_oe)/detC;
            lambda_e = (C_oo - C_oe)/detC;

            COMB.tmpcombined_ccstats=(COMB.add_ccstats_o*lambda_o/sum(1./COMB.ccVars(j_o)) ...
                      + COMB.add_ccstats_e*lambda_e/sum(1./COMB.ccVars(j_e))) ...
                      / (lambda_o + lambda_e) ;

            COMB.tmpcombined_ccspectra = (COMB.add_ccspectra_o*lambda_o/sum(1./ ...
               COMB.ccVars(j_o)) + COMB.add_ccspectra_e*lambda_e/sum(1./COMB.ccVars(j_e))) ...
               / (lambda_o + lambda_e) ;

            COMB.tmpcombined_sensint = (COMB.add_sensint_o + COMB.add_sensint_e - 2*COMB.dataIJ)...
                       *((C_oo*C_ee)/detC);
            COMB.tmperrorBar = sqrt(1/(lambda_o + lambda_e));

          elseif COMB.addctr==2 %only one interval, so no overlapping stuff
            COMB.tmpcombined_ccstats = COMB.add_ccstats_o / sum(1./COMB.ccVars);
            COMB.tmpcombined_ccspectra = COMB.add_ccspectra_o / sum(1./COMB.ccVars);
            COMB.tmpcombined_sensint = COMB.add_sensint_o;
            COMB.tmperrorBar = sqrt(sum(COMB.ccVars))/params.segmentDuration;
          end

          %so, we calculated the combined spectra, ptest etc for the previous
	  %overlapped section of this job, so now add it to 
          %the final combined variables
          if COMB.addctr>=2
            COMB.combined_ccstats = COMB.combined_ccstats + COMB.tmpcombined_ccstats / COMB.tmperrorBar^2;
            COMB.combined_ccspectra = COMB.combined_ccspectra + COMB.tmpcombined_ccspectra / COMB.tmperrorBar^2;
            COMB.combined_sensint = COMB.combined_sensint + COMB.tmpcombined_sensint;
            COMB.ccVars_ovl(COMB.ovladdctr) = COMB.tmperrorBar^2;
            COMB.ovladdctr = COMB.ovladdctr + 1;
          end

          %now, reinitialize the variables for the next overlapping section
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
        end

      else %no overlap

        if params.doBadGPSTimes %determine if the current interval is bad
          if 1
           c1 = params.intervalStartTime - params.bufferSecsMax < COMB.badtimesend;
           c2 = params.intervalStartTime + 3*params.segmentDuration + params.bufferSecsMax > COMB.badtimesstart;
          else
           c1 = params.intervalStartTime+60 == COMB.badtimesstart;
           c2 = 1;
          end

          if sum(c1&c2) > 0 | sqrt(result.ccVar/result.naiVar)>params.maxDSigRatio | sqrt(result.ccVar/result.naiVar)<params.minDSigRatio
            COMB.badtimesflag = 1;
            COMB.badGPSTimes(COMB.badGPSTimesctr) = segmentStartTime;
            COMB.badGPSTimesctr = COMB.badGPSTimesctr + 1;
          else
            COMB.badtimesflag = 0;
          end
        else
          COMB.badtimesflag = 0;
        end

        if ~COMB.badtimesflag
          COMB.add_ccstats = COMB.add_ccstats + result.ccStat / result.ccVar / params.segmentDuration;
          COMB.add_ccspectra = COMB.add_ccspectra + result.ccSpec.data / result.ccVar /params.segmentDuration;
          COMB.add_sensint = COMB.add_sensint + result.sensInt.data * params.segmentDuration^2;
          COMB.ccVars(COMB.addctr) = result.ccVar;
          COMB.addctr = COMB.addctr + 1;
        end
      end %if params.doOverlap
    end %if params.doCombine
