function maptest_rad(plot_file,Trials,matstr)
% function maptest_rad(plot_file)
% uses NEON to calculate max_snr distribution for radiometer map
% plot_file : spherial harmonic plot data file 
% e.g., plot_file = '/archive/home/ethrane/SpH_results/plots.mat'

load(plot_file);

% calculate transforms to go purely real lm-basis ### jan 21 (Stefan and Eric)
[C2R, R2C] = getC2RMatrix(L);

% convert Fisher matrix to real basis ### jan 21 (Stefan and Eric)
fisher0 = C2R * fisher0 * R2C;

% map resolution
res=1;

% dirty sigma map in pixel basis
[d,RA,DECL,dOmg,U]=diagPixel(Fisher,res);

% radiomater sigma map in pixel basis
s_radio = real(d.^-0.5);

% diagonalize Fisher matrix
[V,D] = eig(fisher0);

% eigenvalues of fisher0
egnvls = diag(D);

try, Trials; 
catch, Trials = 100;
end

for ii=1:Trials;
  %fprintf('%i\n', ii);
  % create random vector in natural basis
  z0 = randn(length(D),1);

  % scale by sigmas
  z = sqrt(egnvls).*z0;

  % convert to lm basis
  x = V*z;

  % convert to complex basis ### jan 21 (Stefan and Eric)
  x = R2C*x;

  % dirty map in pixel basis
  [mapr,ra,decl] = makemap(x,res);

  % radiometer map in pixel basis
  p_radio = real(mapr./d); 

  % record max and min snr
  max_snr(ii) = max(p_radio./s_radio); % max radiometer snr
  min_snr(ii) = min(p_radio./s_radio); % minimum radiometer snr

  % location of max radiometer snr
  try
    ra_max(ii) = ra(p_radio./s_radio==max(p_radio./s_radio));
    decl_max(ii) = decl(p_radio./s_radio==max(p_radio./s_radio));
  catch
    % hot spot at one of the poles
    ra_max(ii) = 0;
    decl_max = 100;
  end
end

try, matstr;
catch, matstr='maptest_rad';
end

save([matstr '.mat']);

return;
