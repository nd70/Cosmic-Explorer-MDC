function UL = getUpperLimitFromMC(MC, Y, Y_sig, conf, calErr)
% function UL = getUpperLimitFromMC(MC, Y, Y_sig, conf, calErr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fuction is an extension to getUpperLimit.m,
% i.e. it marginalizes over the calibration error
% and calculates the upper limit.
% But instead of assuming a Gaussian distribution,
% it take a sample data (trials) as first argument.
%
% Input arguments:
%  MC        : array of MC data [Ntrial x n] where n is the number of
%              [Y,sigma] pairs.  E.g., n can run from 0 to lmax
%  Y         : point estimate [n x 1]
%  Y_sig     : error bar [n x 1]
%  conf      : confidence level (e.g. 0.9 for 90%)
%  calErr    : calibration error (e.g. 0.1 for 10%)
%              currently limited to <0.2 (due to integrand singularity)
%
% Output argument
%  UL        : upper limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eric Thrane & Stefan Ballmer 2011 Feb 23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(Y)-1;
for l=0:L
  % Define the domain over which the likelihood function is defined
  Omg = -10*Y_sig(l+1) : Y_sig(l+1)/500 : 10*Y_sig(l+1);

  % testing and debugging flag
  Marginalize = 'true';
%  Marginalize = 'false';

  if Marginalize
    % marginalized likelihood
    LklhdM = zeros(size(Omg));  
%    for lambda=1-3*calErr:calErr/50:1+3*calErr
    for lambda=calErr/50:calErr/50:1+3*calErr
      % likelihood function defined such that the mean is Y
      Lklhd0 = hist(MC(:,l+1)+Y(l+1), Omg);
      Lklhd0=Lklhd0/sum(Lklhd0);

      % likelihood function given lambda
      % ... but now scale mean and sigma by lambda
      Lklhd = hist( (MC(:,l+1)+Y(l+1))/lambda, Omg);
      % kill the overflow
      Lklhd(  1)=0;
      Lklhd(end)=0;

%WRONG      % normalize likelihood function given lambda
%WRONG      Lklhd = Lklhd/sum(Lklhd);

      % define prior-----------------------------------------------------------
      % Gaussian prior (may not formally converge)
      prior=exp(-0.5*(lambda-1)^2/calErr^2);

      % flat prior
%      prior = (heaviside(lambda-1+calErr).*heaviside(calErr-lambda+1)) ...
%        /(2*calErr);
      % heaviside functions produce NaNs at heaviside=0
      % the solution is to make it a gradual step
%      if(isnan(prior)), prior=1/(4*calErr); end;
      %------------------------------------------------------------------------

      % calculate marginalized likelihood
      LklhdM = LklhdM + Lklhd * prior;

      % How to calculate mean and sigma of Lklhd distribution.
      %meanL0 = sum(Omg.*Lklhd0);
      %meanL = sum(Omg.*Lklhd);
      %stdL0 = sqrt( sum(Omg.^2.*Lklhd0) - meanL0^2 );
      %stdL = sqrt( sum(Omg.^2.*Lklhd) - meanL^2 );
    end
  else
    % no marginalization
    % likelihood function defined such that the mean is Y
    Lklhd0 = hist(MC(:,l+1)+Y(l+1), Omg);
    LklhdM = Lklhd0;
  end

  % calculate posterior assuming positivity
  posterior = LklhdM;
  posterior(Omg<0) = 0;

  % normalize
  posterior = posterior/sum(posterior);

  % define cumulative distribution
  cdf = cumsum(posterior,2);

  % find index where cdf crosses conf
  idx = abs(cdf-conf)==min(abs(cdf-conf));

  if sum(idx)>1
    UL(l+1) = sum(Omg(idx))/sum(idx);
    fprintf('test message: %i\n', sum(idx));
  elseif sum(idx)<1
    UL(l+1) = 0;
  else
    UL(l+1) = Omg(idx);
  end
end

return
