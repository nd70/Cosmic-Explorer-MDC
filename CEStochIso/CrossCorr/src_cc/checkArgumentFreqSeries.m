function Hout=checkArgumentFreqSeries(H,flow,deltaF,numFreqs,interpolateLogarithmic)
%  function Hout=checkArgumentFreqSeries(H,flow,deltaF,numFreqs,interpolateLogarithmic)
%
%  checkArgumentFreqSeries  -- converts H into a frequency Series
%                              H can be a filename, Nx2 matrix or a freq Series
%                    
%  arguments: H        - input object
%             flow     - lowest frequency
%             deltaF   - freq. step size
%             numFreqs - number of frequencies
%
%  output:    Hout   - frequency Series
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
  interpolateLogarithmic;
catch
  interpolateLogarithmic = true;
end;

try
  H.flow;
catch
  if size(H,2) ~= 2
    filename=H;
    H=load(filename);
  end;
  f=(0:(numFreqs-1))'*deltaF + flow;
  if interpolateLogarithmic
    if H(1,1)==0
      fprintf('Warning: using logarithmic interpolation - ignoring H=%g at %gHz\n',H(1,2),H(1,1));
      H=H(2:end,:);
    end
    if f(1)==0
      if f(2)<H(1,1)
        error('Need to specify simulated signal shape H down to %gHz (DC is disregarded)',f(2));
      elseif f(end)>H(end,1)
        error('Need to specify simulated signal shape H up to %gHz',f(end));
      end
      data=[0;exp(interp1(log(H(:,1)),log(H(:,2)),log(f(2:end))))];
    else
      if f(1)<H(1,1)
        error('Need to specify simulated signal shape H down to %gHz',f(1));
      elseif f(end)>H(end,1)
        error('Need to specify simulated signal shape H up to %gHz',f(end));
      end
      data=exp(interp1(log(H(:,1)),log(H(:,2)),log(f)));
    end
  else
    if f(1)<H(1,1)
      error('Need to specify simulated signal shape H down to %gHz',f(1));
    elseif f(end)>H(end,1)
      error('Need to specify simulated signal shape H up to %gHz',f(end));
    end
    data=interp1(H(:,1),H(:,2),f);
  end
  if any(isnan(data))
    error('Simulated signal shape H contains NaN');
  end
  H=constructFreqSeries(data, flow, deltaF);
end; % catch

% check that frequency series have the correct length
if ( (length(H.data) ~= numFreqs) )
  error('size mismatch');
end
% check that frequency series have the correct flow
if ( (H.flow ~= flow) )
  error('flow mismatch');
end
% check that frequency series have the same deltaF
if ( (H.deltaF ~= deltaF) )
  error('deltaF mismatch');
end
Hout=H;
