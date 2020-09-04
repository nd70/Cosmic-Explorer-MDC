function printmap(array, xvals, yvals, xlab, ylab, zlab, crange, doFlip, ...
  doLogZ)
% function printmap(array, xvals, yvals, xlab, ylab, zlab, crange, doFlip, ...
%   doLogZ)
% E. Thrane - make a ft map plot, color range is optional
% e.g.: printmap(ftmap_plot,xvals,yvals,'time (s)','f (Hz)','SNR',[-5 5])
% ftmap_plot has indices (y,x)

  if exist('doFlip')
    if doFlip==1
      array=flipdim(array,1); %y-values are flipped
    end
  end
  iptsetpref('ImshowAxesVisible','on');
  figure;
  imshow(array,[],'Xdata',xvals,'Ydata',yvals, ...
	 'InitialMagnification','fit');
  colormap(hot), h=colorbar;
  set(h, 'FontSize', 20);
  daspect([1 1 1]), axis square;
  ylabel(h,zlab)
  xlabel(xlab);
  ylabel(ylab);
  if exist('crange')
    caxis(crange);
  end
  iptsetpref('ImshowAxesVisible','on');
  axis xy;
  axis_val = get(gca,'Clim');
  axis_fact = range(axis_val)/abs(axis_val(1));
  if exist('doLogZ')
    if (doLogZ==1)
      set(h,'Yscale','log');
    end
  else
  end
return
