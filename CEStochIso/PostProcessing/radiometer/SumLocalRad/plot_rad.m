function plot_rad(pproc_params,Y,E,ff,scale_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Patrick Meyers; Modified
% from similar code written by Eric Thrane
% Creates a set of Diagnostic Plots for use
% in Evaluating output of Radiometer Code.
% Currently, this DOES NOT
% calculate upper limits at all. These plots are
% more useful for background/sensitivity studies
% and are produced typically for each
% detector pair and epoch.
% 
% Parameters:
% -----------
%   pproc_params : struct
%       post processing parameters as read in from
%       a paramsFile using readParamsFromFilePostProc
%   Y : float vector
%       Point Estimate / sigma^2
%   E : float vector
%       1/sigma^2
%   ff : float vector
%       frequency vector containing frequencies
%       used in analysis
%   scale factor : float
%       scaling factor from calibration group
%
% CONTACT: patrick.meyers@ligo.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


index = 0;
for mm = 1:(pproc_params.numSkyDirections)
  if any(pproc_params.skippedSkyDirections == mm)
    continue
  end
  index = index+1;
   % beginning of title for plots "H1L1 A" for example

   % define sky direction and name of
   % direction to be used in title for plots
   % associatd with that sky direction
   skydirection  = pproc_params.skyDirectionName{index};
   skydirtitle   = pproc_params.skyDirectionTitle{index};
   save_loc      = pproc_params.output_plot_dir_prefix;
   %str_title     = pproc_params.plotTitlePrefix{index};

   % for the purposes of saving things
   str           = ['_' skydirection '_plots_'];
   % amplitude scaling factor
   path          = [save_loc '_' skydirection];
   % Calculate Y, sigma, and scale them based on scale factor calculated in
   % S6 error budget

   % perform proper scalings now
   % Note that from this point forward
   % Y has units of strain^2
   % and sigma has units of strain
   % because we multiply by deltaF
   %
   % SEE LIGO-T040128-00-E for details on Bias Factor
   N = 2*9/11*(2*pproc_params.segmentDuration*pproc_params.deltaF-1);
   bias_factor = N/(N-1);
   deltaF  = pproc_params.deltaF;
   pte     = E(:,index).^-1.*Y(:,index)*scale_factor*deltaF;
   sig     = E(:,index).^-0.5;
   sig     = deltaF*sig*scale_factor*bias_factor;

   % calculate SNR
   snr = pte./sig;
   idx = find(snr==max(abs(snr))|snr==-max(abs(snr)));
   fprintf('max snr = %2.2f at f=%4.2f\n',snr(idx), ff(idx));
   fprintf('max pte = %1.2e at f=%4.2f\n',sqrt(pte(idx)),ff(idx));


  
	 % cut out extra frequencies
	 % for post processing
   % cut out NaN values
   cut   = ~(isnan(pte)|isnan(sig));
	 
	 
   pte   = pte(cut);
   sig   = sig(cut);
   f     = ff(cut);
   snr   = snr(cut);
	 final_mask = ff(~cut);
   % perform ks test on snr to determine if our distribution looks gaussian
   [h,KSSTAT] = kstest(snr(find(snr~=max(abs(snr)))));
   [MEAN,STD,MEANCI,STDCI] = normfit(snr(find(snr~=max(abs(snr)))));
   MEDIAN = median(snr);
   fprintf('%4.4f is the p-value of the snr distribution for this run\n',KSSTAT);
   fprintf('%4.4f is the standard deviation of the snr distribution\n',STD);
   fprintf('%4.4f is the median snr value\n',MEDIAN);
   fprintf('%4.4f is the mean snr value\n',MEAN);



   % write loud bins (snr>4) to a file
   fid = fopen([save_loc str 'combined_loud_bins.txt'], 'w+');
   fprintf(fid, 'Very loud bins, snr>4.0\n');
   for i=1:length(f)-1
     if abs(snr(i)) > 4.0
       fprintf(fid, '%4.2f %2.3f\n', f(i), snr(i));
     end
   end
   fclose(fid);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Make Plots and Text Files
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if pproc_params.doPlots
      % plot SNR
      figure;
      plot(f, snr);
      xlabel('frequency (Hz)');
      ylabel('snr');
      pretty;
%      axis([48 52 -10 30]);
      title([skydirtitle ': Narrow Band SNR Plot']);
      print('-dpng', [save_loc str 'combined_snr.png']);
      print('-depsc2', [save_loc str 'combined_snr.eps']);

      % plot Y
      figure;
      semilogy(f, abs(pte));
      xlabel('frequency (Hz)');
      ylabel('Point Estimate');
      pretty;
      axis([0 2e4 0 1e-22]);
      title([skydirtitle ': Point Estimate Plot']);
      print('-depsc2', [save_loc str 'combined_ptest.eps']);


      % plot Y and sigma together
      figure;
      loglog(f, abs(pte), 'g');
      hold on;
      loglog(f, sig, 'b');
      xlabel('frequency');
      title([skydirtitle ':Point Est. and  \sigma vs. Frequency ' ]);
      legend('Point Estimate', 'Sensitivity');
      axis([40 1800 0 10^-46 ]);
      pretty;
      print('-depsc2', [save_loc str 'combined_ptest_sig.eps']);
      print('-dpng','-r400', [save_loc str 'combined_ptest_sig']);



      % plot histogram of SNR values. Should be normally distributed, so fit it
      % with a normal distribution to be sure
      figure;
      histfit(snr(find(snr~=max(abs(snr)))), 50, 'normal');
      xlabel('snr');
      ylabel('number of frequency bins');
      title([skydirtitle ': histogram of snr values in frequency bins']);
      pretty;
      print('-dpng', [save_loc str 'combined_snr_hist.png']);
      print('-depsc2', [save_loc str 'combined_snr_hist.eps']);

      figure;
      semilogy(f,sqrt(sig),'b');
      xlabel('frequency (Hz)');
      ylabel('strain in each bin');
      title([skydirtitle ': strain in each bin']);
     % axis([48 52 4*10^-25 10^-22]);
      grid on;
      pretty;
      print('-depsc2',[save_loc str 'combined_sig.eps']);

   end

   %%%%%%%%%%%%%%%%%%%%%%%
   % Save Final Data
   %%%%%%%%%%%%%%%%%%%%%%%
   % This data is scaled and Y is an estimate of total
   % strain in each bin (units of strain^2)
   % and sigma is the strain noise in each bin 
   % (units ofstrain^2) as well
   %%%%%%%%%%%%%%%%%%%%%%%

   ptEst = pte;


   save([path '_final_combined_data.mat'],'ptEst','sig','f','snr', 'STD', 'KSSTAT','final_mask');
end
