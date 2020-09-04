function [asigthresh]=set_absolute_sigma_cut(ifopair,flow,epoch);
%keyboard
% sets absolute sigma cut based on the interferometer pair and frequency range.
% the cuts decided upon below are supp|ted in the S6 radiometer documentation.
% ifopair = H1L1, H1V1, L1V1
% flow = 41.5 | flow == 40, 170, 600, 1000
if strcmp(ifopair,'H1L1')
  if strcmp(epoch,'A')
    if flow == 41.5 | flow == 40
      asigthresh = 2;
    elseif flow == 170
      asigthresh = 1.25;
    elseif flow == 600
      asigthresh = 3.75;
    elseif flow == 1000
      asigthresh = 8;
    end
  elseif strcmp(epoch,'B')
    if flow == 41.5 | flow == 40 
      asigthresh = 2;
    elseif flow == 170
      asigthresh = 1;
    elseif flow == 600
      asigthresh = 4.25;
    elseif flow == 1000
      asigthresh = 8;
    end
  elseif strcmp(epoch,'C')
    if flow == 41.5 | flow == 40
      asigthresh = 1.75;
    elseif flow == 170
     asigthresh = 0.8;
    elseif flow == 600
      asigthresh = 3;
    elseif flow == 1000
      asigthresh = 6;
    end
  elseif strcmp(epoch,'D')
    if flow == 41.5 | flow == 40
      asigthresh = 2;
    elseif flow == 170
     asigthresh = 0.6;
    elseif flow == 600
      asigthresh = 2.5;
    elseif flow == 1000
      asigthresh = 5;
    end
  end
elseif strcmp(ifopair,'H1V1')
  if strcmp(epoch,'A')
    if flow == 41.5 | flow == 40
      asigthresh = 4;
    elseif flow == 170
      asigthresh = 2;
    elseif flow == 600
      asigthresh = 3.5;
    elseif flow == 1000
      asigthresh = 6;
    end
  elseif strcmp(epoch,'B')
    if flow == 41.5 | flow == 40
      asigthresh = 5;
    elseif flow == 170
      asigthresh = 2;
    elseif flow == 600
      asigthresh = 6;
    elseif flow == 1000
      asigthresh = 11;
    end
  elseif strcmp(epoch,'D')
    if flow == 41.5 | flow == 40
      asigthresh = 5;
    elseif flow == 170
     asigthresh = 2;
    elseif flow == 600
      asigthresh = 6;
    elseif flow == 1000
      asigthresh = 9;
    end
  end
elseif strcmp(ifopair,'L1V1')
  if strcmp(epoch,'A')
    if flow == 41.5 | flow == 40
      asigthresh = 5;
    elseif flow == 170
      asigthresh = 2.25;
    elseif flow == 600
      asigthresh = 5;
    elseif flow == 1000
      asigthresh = 8.5;
    end
  elseif strcmp(epoch,'B')
    if flow == 41.5 | flow == 40
      asigthresh = 6;
    elseif flow == 170
      asigthresh = 3.5;
    elseif flow == 600
      asigthresh = 8;
    elseif flow == 1000
      asigthresh = 12;
    end
  elseif strcmp(epoch,'D')
    if flow == 41.5 | flow == 40
      asigthresh = 6;
    elseif flow == 170
     asigthresh = 2;
    elseif flow == 600
      asigthresh = 8;
    elseif flow == 1000
      asigthresh = 11;
    end
  end
end
asigthresh = asigthresh*10^-45;
