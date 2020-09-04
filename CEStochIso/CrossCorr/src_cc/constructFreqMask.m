function mask = constructFreqMask(fMin, fMax, deltaF, ...
                                  freqsToRemove, nBinsToRemove, doFreqMask)
%
%  constructFreqMask --- constructs frequency mask 
%
%  constructFreqMask(fMin, fMax, deltaF, freqsToRemove, nBinsToRemove, 
%  doFreqMask) returns an array of 1's and 0's used to ignore certain
%  frequencies in the cross-correlation sum.
% 
%  freqsToRemove is an array of discrete frequency bins to remove.
%
%  nBinsToRemove is an array containing the number of frequency bins to 
%  remove for the corresponding discrete frequency:
%
%  -  <=0 means don't remove that particular line
%  -  if odd, do symmetric removal of bins
%  -  if even, remove one more high frequency bin than low frequency bin
%         
%  doFreqMask = true/false to construct/not construct frequency mask  
%
%  Routine written by Joseph D. Romano
%  Contact Joseph.Romano@astro.cf.ac.uk
%
% $Id: constructFreqMask.m,v 1.4 2005-02-24 13:03:39 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFreqs = floor((fMax-fMin)/deltaF)+1;
mask = ones(numFreqs,1);  % trivial case

% check that there are frequencies to remove
if isempty(freqsToRemove) 
  return 
end

% construct frequency mask
if doFreqMask

%%%%%%% Commented out, sGc
  % check that frequency range is okay ...
%  if ( (freqsToRemove(1) < fMin) | (freqsToRemove(end) > fMax) )
%    error('invalid frequencies to remove');
%  end
%%%%%%%

  % excise bad frequencies
  for ii=1:length(freqsToRemove)

    if nBinsToRemove(ii)>0    

      index = floor( (freqsToRemove(ii)-fMin)/deltaF ) + 1;

      if mod(nBinsToRemove(ii),2)==0
        % even number of bins to remove
        index_low  = index - nBinsToRemove(ii)/2 + 1;
        index_high = index + nBinsToRemove(ii)/2;
      else 
        % odd number of bins to remove
        index_low  = index - (nBinsToRemove(ii)-1)/2;
        index_high = index + (nBinsToRemove(ii)-1)/2;
      end

      % truncate index values if they are too small or too large
      if index_low < 1 
        index_low=1;
      end

      if index_high > numFreqs 
        index_high = numFreqs;
      end

      % add zeros to frequency mask
      if index_high >= index_low % added, sGc
        mask(index_low:index_high) = zeros(index_high-index_low+1,1);
      end % added, sGc

      % print some warnings if necessary
      if (freqsToRemove(ii) < fMin | freqsToRemove(ii) > fMax)
        fprintf('WARNING: Requested to notch frequency (%g Hz) not in frequency band (%g - %g Hz)\n',freqsToRemove(ii),fMin,fMax)
      end

    end % if nBinsToRemove(ii)>0

  end % loop over freqs to remove

end

return

