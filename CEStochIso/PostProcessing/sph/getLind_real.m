function [ind] = getLind_real(Lmax, l)
  
  if l==0
    ind = 1;
  else
    ind = [];
    for m=0:l
      for k=0:1
        ind = [ind ; 1 + l^2 + max(2*m-1,0) + k];
      end
    end
    ind=unique(ind);
  end
return
