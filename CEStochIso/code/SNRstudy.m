load fullComb_diffFreqMask_ptEstIntegrand.dat
load fullComb_diffFreqMask_sensIntegrand.dat
ptEstIntegrand=fullComb_diffFreqMask_ptEstIntegrand(:,3);
sensIntegrand=fullComb_diffFreqMask_sensIntegrand(:,3);
sens_f=sensIntegrand/32;
sens_total=sum(sens_f);
ptEst_f=2/32*sens_total*ptEstIntegrand./sens_f;
freq=fullComb_diffFreqMask_sensIntegrand(:,2);




sigma_f=sqrt(sens_f.^-1)*1.06;
nocut=isfinite(sigma_f);
snr_f=ptEst_f(nocut)./sigma_f(nocut);
f=freq(nocut);

cum_sens_total=cumsum(sens_f(nocut));
cum_sigma=sqrt(cum_sens_total.^-1)*1.06;
cum_ptEst=cumsum(sens_f(nocut).*ptEst_f(nocut))./cum_sens_total;
cum_snr_f=cum_ptEst./cum_sigma;

h=figure;
addpath(genpath('/home/sgwynne.crowder/bin/matlab/'))
pretty
plot(f,snr_f,'k')
hold on
plot(f,cum_snr_f,'r')
xlabel('Frequency (Hz)')
ylabel('SNR')
axis([20 100 -5 5])
title('192s a0, red = cum SNR')
%print(h,'-dpng','cumSNR_192s_a0_fullCombNotch.png')

cumsens=cumsum(sensIntegrand(nocut));
cumsens=cumsens/cumsens(end);

histfit(snr_f)
fitdist(snr_f,'normal')
