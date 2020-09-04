function errorsvsLmax(detectorPairs, Lstop, flow, fhigh, alpha)
%
% calculate errors in PLM as a function of lmax for a set of detector pairs
%
% Input:
% 
%   detectorPairs - a cell array containing a set of strings
%                   corresponding to the following detector pairs 
%
%        HL: Hanford-Livingston
%        HV: Hanford-Virgo
%        LV: Livingston-Virgo
%
%   Lstop  - final value of L (e.g., 15)
%   flow   - low frequency (Hz) (e.g., 50 Hz)
%   fhigh  - high frequency (Hz) (e.g., 1000 Hz)
%   alpha  - exponent in stochastic background power law (e.g, 0 or -3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
warning('off','MATLAB:log:logOfZero');
close all
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontWeight','bold');
set(0,'DefaultTextFontSize',16);
%set(0,'DefaultTextFontWeight','normal');
%set(0,'DefaultUicontrolFontSize',12);
%set(0,'DefaultUicontrolFontWeight','normal');
%set(0,'DefaultAxesTickLength',[0.02 0.01]);
set(0,'DefaultAxesLineWidth', 2.0);
set(0,'DefaultLineLineWidth', 2.0);
%set(0,'DefaultAxesGridLineStyle','--');
set(0,'DefaultAxesXColor',[0 0 0]);
set(0,'DefaultAxesYColor',[0 0 0]);

numBaselines = numel(detectorPairs);

% frequency array
numFreqs = 7000;
deltaF = (fhigh-flow)/numFreqs;
f = transpose(flow + deltaF*[0:numFreqs-1]);

% initialise variables
P1 = zeros(numFreqs,numBaselines);
P2 = zeros(numFreqs,numBaselines);
H = zeros(numFreqs);

for ii=1:numBaselines

  % get SRD data (f, amplitude spectral density) for the various detectors
  switch detectorPairs{ii}
    case 'HL'
      %srddata1=load('LIGOsrd.dat');
      %srddata1(:,2)=srddata1(:,2)/4e3; % convert LIGO srd to strain
      %srddata2=load('LIGOsrd.dat');
      %srddata2(:,2)=srddata2(:,2)/4e3; % convert LIGO srd to strain
      srddata1=load('lho4k_070318_strain.txt');
      srddata2=load('llo_060604_strain.txt');
    case 'HV'
      %srddata1=load('LIGOsrd.dat');
      %srddata1(:,2)=srddata1(:,2)/4e3; % convert LIGO srd to strain
      %srddata2=load('VIRGOsrd.dat');
      srddata1=load('lho4k_070318_strain.txt');
      srddata2=load('SensitivityH_VSR1.txt');
    case 'LV'
      %srddata1=load('LIGOsrd.dat');
      %srddata1(:,2)=srddata1(:,2)/4e3; % convert LIGO srd to strain
      %srddata2=load('VIRGOsrd.dat');
      srddata1=load('llo_060604_strain.txt');
      srddata2=load('SensitivityH_VSR1.txt');
    otherwise
      error('not a valid detector pair');
  end

  % calculate power spectra, H(f)
  fprintf('Calculating PSDs for detector pair %s\n', detectorPairs{ii});
  srd1=interp1(srddata1(:,1),srddata1(:,2),f);
  srd2=interp1(srddata2(:,1),srddata2(:,2),f);
  P1(:,ii)=srd1.^2;
  P2(:,ii)=srd2.^2;
  H = f.^alpha;

end

% initialise variables
errors = cell(Lstop+1,Lstop+1);

% calculate covariance matrix and error in PLMs for different Lmax
for Lmax=0:Lstop

  % initialise variables
  fisher_matrices = zeros((Lmax+1)^2,(Lmax+1)^2,numBaselines);

  fprintf('\n');
  fprintf('Working on Lmax=%d out of %d\n',Lmax,Lstop); 

  for ii=1:numBaselines

    % calculate gammaLMs
    fprintf('Calculating gammaLMs for detector pair %s out to Lmax=%d\n', ...
            detectorPairs{ii}, Lmax);
    glm = calGammaLM('', detectorPairs{ii}, Lmax, numFreqs, flow, deltaF);
    lvec = glm.lvec;
    mvec = glm.mvec;

    % construct fisher information matrix
    fprintf('Constructing Fisher matrix for detector pair %s\n', ...
            detectorPairs{ii});
    mask = ones(length(P1(:,ii)),1);
    fisher_matrices(:,:,ii)=calFisher(P1(:,ii), P2(:,ii), H, glm, 0, 60, 1, 1, mask, 1);

  end;

  % construct fisher matrix for multiple baselines
  fisher_matrix = zeros((Lmax+1)^2,(Lmax+1)^2);
  for ii=1:numBaselines
    fisher_matrix = fisher_matrix + fisher_matrices(:,:,ii);
  end 

  % invert fisher matrix to get covariance matrix 
  % (note: use invLMblock to ignore the non-zero m!=m' elements in the
  % fisher matrix.  these elements will equal zero when integrating over
  % one sidereal day for time-independent power spectra.)
  fprintf('Inverting Fisher matrix\n');
  cov_matrix = invLMblock(fisher_matrix);
  %cond(cov_matrix)

  % make image plot of log(abs(covariance)) matrix
  if 1
    figure(Lmax+1);
    imagesc(log(abs(cov_matrix)));
    colorbar
    fig_title=['Log(abs(cov matrix)) for Lmax = ' num2str(Lmax)];
    title(fig_title)

    fname=['covMatrix' num2str(Lmax) '.jpg'];
    print('-djpeg',fname);

  end

  % calculate errors in PLM estimates for M=0,1,...Lmax
  for M=0:Lmax  
    ind=find(mvec==M);
    ii = 0;
    for L=lvec(ind(1)):lvec(ind(end))
      % diagonal elements are real and positive (to machine round-off)
      errors{L+1,M+1} = [errors{L+1,M+1}, ...
        sqrt(real(diag(cov_matrix(ind(1)+ii,ind(1)+ii))))];
      ii = ii+1;
    end
  end

end

% plot error in P00, P10, P20, P30 as a function of Lmax
figure(Lstop+2);
semilogy(0:Lstop,errors{1,1},'*-b')
hold on
semilogy(1:Lstop,errors{2,1},'*-g')
hold on
semilogy(2:Lstop,errors{3,1},'*-r')
hold on
semilogy(3:Lstop,errors{4,1},'*-c')
xlabel('Lmax');
ylabel('Error in plm estimates');
grid on
ylim([1e-46 1e-44])
legend('monopole','dipole','quadrupole', 'octupole','Location','NorthWest');
print('-depsc2','errorPlmvsLmax');

return

