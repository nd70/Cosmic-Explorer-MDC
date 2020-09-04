function data=periodicSubtract(data,correction)
%  function data=periodicSubtract(data,correction)
%
%  Subtracts a periodic version of correction from data
%  This routine is used to subtract the average timing transient
%  in H1 data.
%
%  data       = data input and data output
%  correction = correction to be subtracted from data
%               length(data) must be a multiple of length(correction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dl=length(data);
cl=length(correction);
if mod(dl,cl)==0
  N=dl/cl;
  ind=1:cl;
  for ii=0:N-1
    data(ii*cl+ind)=data(ii*cl+ind)-correction;
  end
else
 fprintf('Data size: %d\n',dl);
 fprintf('Corr size: %d\n',cl);
 error('data must be a multple of correction');
end
return
