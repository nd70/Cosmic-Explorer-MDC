function maptest_clean(plot_file,Trials,matstr)
% function maptest(plot_file)
% e.g., plot_file = '/archive/home/ethrane/SpH_results/plots.mat';
% uses NEON to calculate max_snr distribution for radiometer map

% load covariance matrix mat file
load(plot_file);

% calculate transforms to go purely real lm-basis ### jan 21 (Stefan and Eric)
  [C2R, R2C] = getC2RMatrix(L);

% convert pCovar0 matrix to real basis ### jan 21 (Stefan and Eric)
%  fisher0 = C2R * fisher0 * R2C;
  pCovar0 = C2R * pCovar0 * R2C;

% map resolution
res=1;

% pCovar0 sigmaPix0 map are already calculated by plots.m

% diagonalize covariance matrix
[V,D] = eig(pCovar0);

% eigenvalues of Fisher
egnvls = diag(D);

try, Trials; 
catch, Trials = 1000;
end

for ii=1:Trials
  %fprintf('%i\n', ii);
  % create random vector in REAL natural basis
  z0 = randn(length(pCovar0),1);

  % scale by sigmas
  z = sqrt(egnvls).*z0;

  % convert to lm basis
  p = V*z;

  % convert to complex basis ### jan 21 (Stefan and Eric)
  p = R2C*p;

  % clean map in pixel basis
  [mapr,ra,decl] = makemap(p,res);

  % record max and min snr
  max_snr(ii) = max(mapr./sigmaPix0); % max clean map snr
  min_snr(ii) = min(mapr./sigmaPix0); % minimum clean map snr

  % location of max radiometer snr
  try
    ra_max(ii) = ra(mapr./sigmaPix0==max(mapr./sigmaPix0));
    decl_max(ii) = decl(mapr./sigmaPix0==max(mapr./sigmaPix0));
  catch
    % hot spot at one of the poles
    ra_max(ii) = 0;
    decl_max = 100;
  end
end

try, matstr;
catch, matstr='maptest_clean';
end

save([matstr '.mat']);

return;
