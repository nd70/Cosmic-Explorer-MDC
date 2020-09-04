function [x_tot, fisher_tot, nbad] = combineXFisher(vSpH, sig, naiSig, params)
  % by E. Thrane on Sept. 1, 2010
  % initialize variables
  Q = (vSpH.Lmax+1)^2;
  x_tot=zeros(Q,1);  fisher_tot=zeros(Q,Q);
  fisher1=zeros(Q,Q);  fisherN=zeros(Q,Q);
  epsilon = 3/70;   % factor from overlapping Hann windows
  nbad = 0; % number of segments failing sigma cut

  % the zeroth segment is always considered bad in order to set fisher1
  lastsegmentwasbad=1;

  for ii=1:length(vSpH.out.data)
%fprintf('%i:\n', ii);
    if sig(ii)/naiSig(ii)<params.minDSigRatio | ...      % This segment is bad.
         sig(ii)/naiSig(ii)>params.maxDSigRatio
      nbad = nbad + 1;
      lastsegmentwasbad=1;
%fprintf(' segment %i is bad\n', ii);
    else                                                % This segment is good.
      x_tot = x_tot + vSpH.out.data{ii}.X;
      fisher_tot = fisher_tot + vSpH.out.data{ii}.Fisher;
      % The last segment is bad but this one is good. -------------------------
      if lastsegmentwasbad
%fprintf(' last segment was bad, but this one is good.\n');
%fprintf(' Fisher1=%i and FisherN=%i\n',fisher1,fisherN);
        fisher_tot = fisher_tot+((epsilon/2)/(1-2*epsilon))*(fisher1+fisherN);
        fisher1 = vSpH.out.data{ii}.Fisher;
      end % -------------------------------------------------------------------
      fisherN = vSpH.out.data{ii}.Fisher;
      lastsegmentwasbad=0;
    end
  end
%fprintf('# Apply final corrections:');
%fprintf(' Fisher1=%i and FisherN=%i\n',fisher1,fisherN);
  % apply final epsilon corrections
  x_tot = (1-2*epsilon)*x_tot;
  fisher_tot = (1-2*epsilon)*fisher_tot + (epsilon/2)*(fisher1+fisherN);

return
