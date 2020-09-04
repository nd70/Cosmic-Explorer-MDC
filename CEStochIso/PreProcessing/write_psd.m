
x=load('ZERO_DET_high_P.txt');

fid=fopen('ALIGO-HighP_PSD.txt','w+');

for i=1:length(x(:,1))
   fprintf(fid,'%0.7e %0.7e\n',x(i,1),x(i,1)^2)
end

