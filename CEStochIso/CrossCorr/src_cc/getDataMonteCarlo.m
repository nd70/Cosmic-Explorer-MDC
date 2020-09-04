function [data1,data2]=getDataMonteCarlo(params,data1,data2,I,J,K)

% simulate stochastic signals for the duration of the whole job
% Previously done in the old stochastic.m
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: getDataMonteCarlo.m,v 1.1 2007-07-18 01:18:20 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global MCARLO;


        % determine seed for MC time-offset (indexed by segment number)
        if params.doOverlap
          segmentNumber = (I-1)+(2*J-1);
        else
          segmentNumber = (I-1)+J; 
        end
        seed = params.jobNumber*1000000+(K-1)*params.numSegmentsTotal+segmentNumber;

        %fprintf('job = %d, I = %d, J = %d, K = %d, seg number = %d, seed = %d\n',params.jobNumber,I,J,K,segmentNumber,seed);

        if params.doMCoffset % time-offset detector1
	  if params.doConstTimeShift
	    MCoffset = round(params.ConstTimeShift / data1.deltaT);
          else
	    rand('state',20*seed);
            MCoffset = round(((params.maxMCoff-params.minMCoff) * rand(1) + params.minMCoff) /...
                             data1.deltaT);
          end

          tempdata1 = circshift(data1.data,MCoffset);
        else
          tempdata1 = data1.data;
        end % if params.doMCoffset

        % inject simulated signal into the detector noise
        startIndex1 = 1 + ...
          (data1.tlow - params.centeredStartTime + params.bufferSecs1)*params.resampleRate1;
        endIndex1 = startIndex1 + ...
          (params.segmentDuration + 2*params.bufferSecs1)*params.resampleRate1 - 1;
        startIndex2 = 1 + ...
          (data2.tlow - params.centeredStartTime + params.bufferSecs2)*params.resampleRate2;
        endIndex2 = startIndex2 + ...
          (params.segmentDuration + 2*params.bufferSecs2)*params.resampleRate2 - 1;

        %fprintf('centered start time = %d, data tlow = %d, start index = %d, end index = %d\n', params.centeredStartTime, data1.tlow, startIndex1, endIndex1);

        % note that detector 1 was shifted, hence uses tempdata1 instead of data1.data
        data1.data= tempdata1 + sqrt(params.simOmegaRef1)*transpose(MCARLO.h1(K, startIndex1:endIndex1));
        data2.data=data2.data + sqrt(params.simOmegaRef2)*transpose(MCARLO.h2(K, startIndex2:endIndex2));


        % for debugging --->
        %if K==1
        %  times1 = params.centeredStartTime - params.bufferSecs1 + [0:Nsamples1-1];
        %  times2 = params.centeredStartTime - params.bufferSecs2 + [0:Nsamples2-1];
        %  inject1 = [zeros(1,startIndex1-1) ...
        %             MCARLO.h1(K,startIndex1:endIndex1) ...
        %             zeros(1,Nsamples1 - endIndex1)];
        %  inject2 = [zeros(1,startIndex2-1) ...
        %             MCARLO.h2(K,startIndex2:endIndex2) ...
        %             zeros(1,Nsamples2 - endIndex2)];
        %               
        %  figure(1);
        %  subplot(2,1,1); plot(times1, MCARLO.h1(K,:), 'b', times1, inject1, 'r');
        %  title('Time series','FontSize',12)
        %  ylabel('Amplitude 1','FontSize',12);
        %                                                                       
        %  subplot(2,1,2); plot(times2, MCARLO.h2(K,:), 'b', times2, inject2, 'r');
        %  xlabel('time (sec)','FontSize',12);
        %  ylabel('Amplitude 2','FontSize',12);
        %end
        % for debugging <---
