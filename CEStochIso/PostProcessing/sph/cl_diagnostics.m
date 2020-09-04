function cl_diagnostics(C, Cl, Cl_sig, Cl_UL)
% function cl_diagnostics(C, Cl, Cl_sig)
% make plots for Cl limit diagnostics

L = length(Cl)-1;
for l=0:L
  % q must go up to at least 10*sigma to match getUpperLimit
  if Cl_sig(l+1)>0,
    q = -10*Cl_sig(l+1) : Cl_sig(l+1)/100 : 10*Cl_sig(l+1); 
  else
    q=linspace(-10*min(Cl_sig),10*max(Cl_sig),100);     
  end

  % calculate likelihood function
  Lklhd = hist(C(:,l+1)+Cl(l+1), q);
  Lklhd=Lklhd/sum(Lklhd);

  % no-signal PDF just to make plots
  pdf = hist(C(:,l+1), q);
  pdf=pdf/sum(pdf);

  % diagnostic plot 1
  if l==5
    figure;
    plot(q, Lklhd, 'r');
    hold on;
    plot(q, pdf, 'b');
    hold on;
    plot([Cl(l+1) Cl(l+1)], 1.2*[0 max(pdf)], 'g');
    hold on;
    plot([Cl(l+1) Cl(l+1)]-Cl_sig(l+1), 1.2*[0 max(pdf)], 'k');
    hold on;
    plot([Cl(l+1) Cl(l+1)]+Cl_sig(l+1), 1.2*[0 max(pdf)], 'k');
    hold on;
    plot([Cl_UL(l+1) Cl_UL(l+1)], 1.2*[0 max(pdf)], 'k--');
    legend('likelihood', 'no signal pdf', 'Cl', '-\sigma', '+\sigma','UL');
    axis([-1.2e-97 1.2e-97 0 8e-3]);
    pretty;
    print('-djpeg', 'cl_test1.jpg');
  end
end

% diagnostic plot 2
figure;
errorbar([0:L], Cl, Cl_sig,'bx')
hold on;
plot([0:L], Cl_UL, 'g');
pretty;
legend('data', 'UL')
%print('-djpeg', 'cl_test2.jpg');

