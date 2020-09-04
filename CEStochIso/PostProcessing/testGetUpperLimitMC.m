
Y=(-2:0.2:2)';
sig=ones(21,1);

data =randn(1e5,1);
dataB=randn(1e5,21);

muOrig     =getUpperLimit  (        Y,sig,0.9,0.19);
muNewAnalyt=getUpperLimitMC('Gauss',Y,sig,0.9,0.19);
muNewMC1e5 =getUpperLimitMC(data   ,Y,sig,0.9,0.19);
muNewMC1e5B=getUpperLimitMC(dataB  ,Y,sig,0.9,0.19);

smartfigure('1 ORDER 1:2');
%smartfigure('POSITION [10,10,600,400]');
smartfigure;
plot(Y,muOrig*0,Y,muNewAnalyt,Y,muNewMC1e5,Y,muNewMC1e5B);
grid on;
xlabel('Y');
ylabel('UL');
legend('muOrig','muNewAnalyt','muNewMC1e5','muNewMC1e5B');

