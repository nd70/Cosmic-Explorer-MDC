function MCARLO=initMonteCarlo(params)

% get the simulated stochastic data
% Previously done in the old stochastic.m
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: initMonteCarlo.m,v 1.1 2007-07-18 01:18:20 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MCARLO;

  % simulate stochastic signals for the duration of the whole job

  % number of samples (including initial and final buffer seconds)
  MCARLO.Nsamples1 = (params.M*params.segmentDuration+2*params.bufferSecs1)*params.resampleRate1;
  MCARLO.Nsamples2 = (params.M*params.segmentDuration+2*params.bufferSecs2)*params.resampleRate2;

  % initialize arrays to hold simulated SB signals
  MCARLO.h1 = zeros(params.numTrials, MCARLO.Nsamples1);
  MCARLO.h2 = zeros(params.numTrials, MCARLO.Nsamples2);

  % loop over number of trials
  for K=1:params.numTrials

    fprintf('Simulating data for trial %d\n', K);

    % put zeros in the initial and final buffer seconds
    Nparams.bufferSecs1 = params.bufferSecs1*params.resampleRate1;
    Nparams.bufferSecs2 = params.bufferSecs2*params.resampleRate2;

    MCARLO.h1(K, 1:Nparams.bufferSecs1) = zeros(1,Nparams.bufferSecs1);
    MCARLO.h2(K, 1:Nparams.bufferSecs2) = zeros(1,Nparams.bufferSecs2);

    MCARLO.h1(K, end-Nparams.bufferSecs1+1:end) = zeros(1,Nparams.bufferSecs1);
    MCARLO.h2(K, end-Nparams.bufferSecs2+1:end) = zeros(1,Nparams.bufferSecs2);

    % loop over number of segments in the job
    for L=1:params.M

      % calculate transfer functions appropriate for this and half of the 
      % next segment (need both of these response functions for splicing)

      if ( strncmp(params.alphaBetaFile1,   'none', length(params.alphaBetaFile1))   | ...
           strncmp(params.calCavGainFile1,  'none', length(params.calCavGainFile1))  | ...
           strncmp(params.calResponseFile1, 'none', length(params.calResponseFile1)) )

        % the data is already calibrated
        transfer1a = constructFreqSeries(ones(params.numFreqs,1), params.flow, params.deltaF);
        transfer1b = transfer1a;

      else
        % get response functions for the appropriate epochs
        calibsec1a = params.centeredStartTime + (L-1)*params.segmentDuration;
        calibsec1b = calibsec1a + floor(params.segmentDuration/2);

        [R1a, responseOK1a] = ...
          calculateResponse(params.cal1.t, params.cal1.f, params.cal1.R0, params.cal1.C0,...
	   params.cal1.alpha, params.cal1.gamma, calibsec1a, params.ASQchannel1);

        [R1b, responseOK1b] = ...
          calculateResponse(params.cal1.t, params.cal1.f, params.cal1.R0, params.cal1.C0,...
	   params.cal1.alpha, params.cal1.gamma, calibsec1b, params.ASQchannel1);
 
        % if either response function is bad, set flag and exit loop
        if (responseOK1a == false) | (responseOK1b == false)
          MCARLO.badResponse = true;
          fprintf('bad response for detector 1 in MC simulation\n');
          continue; % proceed to simulating data for the next segment
        end
  
        % convert to transfer functions (units: counts/strain) 
        transfer1a = convertResponse(params.cal1.f, R1a, params.flow, params.deltaF, params.numFreqs, 1, 0);
        transfer1b = convertResponse(params.cal1.f, R1b, params.flow, params.deltaF, params.numFreqs, 1, 0);
      end

      if ( strncmp(params.alphaBetaFile2,   'none', length(params.alphaBetaFile2))   | ...
           strncmp(params.calCavGainFile2,  'none', length(params.calCavGainFile2))  | ...
           strncmp(params.calResponseFile2, 'none', length(params.calResponseFile2)) )

        % the data is already calibrated
        transfer2a = constructFreqSeries(ones(params.numFreqs,1), params.flow, params.deltaF);
        transfer2b = transfer2a;

      else
        % get response functions for the appropriate epochs
        calibsec2a = params.centeredStartTime + (L-1)*params.segmentDuration;
        calibsec2b = calibsec2a + floor(params.segmentDuration/2);

        [R2a, responseOK2a] = ...
          calculateResponse(params.cal2.t, params.cal2.f, params.cal2.R0, params.cal2.C0,...
	   params.cal2.alpha, params.cal2.gamma, calibsec2a, params.ASQchannel2);

        [R2b, responseOK2b] = ...
          calculateResponse(params.cal2.t, params.cal2.f, params.cal2.R0, params.cal2.C0,...
	   params.cal2.alpha, params.cal2.gamma, calibsec2b, params.ASQchannel2);

        % if either response function is bad, set flag and exit loop
        if (responseOK2a == false) | (responseOK2b == false)
          MCARLO.badResponse = true;
          fprintf('bad response for detector 2 in MC simulation\n');
          continue; % proceed to simulating data for the next segment
        end
  
        % convert to transfer functions (units: counts/strain) 
        transfer2a = convertResponse(params.cal2.f, R2a, params.flow, params.deltaF, params.numFreqs, 1, 0);
        transfer2b = convertResponse(params.cal2.f, R2b, params.flow, params.deltaF, params.numFreqs, 1, 0);
      end

      % simulate the stochastic signals 
	N1 = params.segmentDuration*params.resampleRate1;
        N2 = params.segmentDuration*params.resampleRate2;
        deltaT1 = 1/params.resampleRate1;
        deltaT2 = 1/params.resampleRate2;

        if isfield(params, 'powerIndex') 
          singalType = params.powerIndex;
        elseif isfield(params, 'signalType')
          singalType = params.signalType;
        else
          error('Signal spectrum type not provided \n');
        end

        [h1a, h2a] = simulateSB(0, deltaT1, deltaT2, ...
                              N1, N2,...
                              singalType, params.detector1, params.detector2, ...
			      transfer1a, transfer2a, ...
                              params.nResample1, params.betaParam1, ...
                              params.nResample2, params.betaParam2, ...
                              params.fbase1, params.fbase2);
        [h1b, h2b] = simulateSB(0, deltaT1, deltaT2, ...
                              N1, N2,...
                              singalType, params.detector1, params.detector2, ...
			      transfer1b, transfer2b, ...
                              params.nResample1, params.betaParam1, ...
                              params.nResample2, params.betaParam2, ...
                              params.fbase1, params.fbase2);

    
      % splice the data
      data1 = spliceData(1, transpose(h1a.data), transpose(h1b.data), N1);
      data2 = spliceData2(2, transpose(h2a.data), transpose(h2b.data), N2);
                                                                                
      % fill in the appropriate part of the signal array
      MCARLO.h1(K, 1+(L-1)*N1 + Nparams.bufferSecs1 : L*N1 + Nparams.bufferSecs1) = data1;
      MCARLO.h2(K, 1+(L-1)*N2 + Nparams.bufferSecs2 : L*N2 + Nparams.bufferSecs2) = data2;
                                                                               
    end % loop over segments L

    % for debugging --->
    %figure(1);
    %subplot(2,1,1); plot(MCARLO.h1(K,:));
    %title('Time series','FontSize',12)
    %ylabel('Amplitude 1','FontSize',12);
                                                                                
    %subplot(2,1,2); plot(MCARLO.h2(K,:));
    %xlabel('time (sec)','FontSize',12);
    %ylabel('Amplitude 2','FontSize',12);
    % for debugging <---
                                                                                
  end % loop over trials K

  fprintf('Done with MC simulations\n');
    
