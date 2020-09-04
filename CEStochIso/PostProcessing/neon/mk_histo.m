function [p,q,n] = mk_histo(str, thresh, plotstr)
% function (str, thresh)

th = strassign(thresh);
load(str);

try, plotstr;
catch, plotstr='histo';
end

% bin data
q = 0.1:0.10:6;
n = hist(max_snr, q);

% calculate p
p = sum(max_snr>=th)/length(max_snr);

% plot results
 figure;
 errorbar(q, n, sqrt(n), 'x-');
 hold on;
 plot([th th], [0 max(n)]);
 xlabel('max snr');
 ylabel('counts');
 tit = ['p(max snr>=' num2str(th) ') = ' num2str(100*p) '%'];
 title(tit);
% pretty;
% print('-djpeg', [plotstr '.jpg']);
% print('-depsc2', [plotstr '.eps']);
try, print('-dpng', [plotstr '.png']);, catch, end


%{
q = -6:0.10:-1;
n = hist(min_snr, q);
% plot results
figure;
errorbar(q, n, sqrt(n), 'bo-');
hold on;
plot([-th -th], [0 max(n)], 'r');
xlabel('min snr');
ylabel('counts');
p = sum(min_snr<=-th)/length(min_snr);
tit = ['p(min snr>=' num2str(th) ') = ' num2str(100*p) '% (' num2str(Trials)...
  ' trials)'];
title(tit);
pretty;
print('-djpeg', 'histo_min.jpg');
print('-depsc2', 'histo_min.eps');
 %}
return
