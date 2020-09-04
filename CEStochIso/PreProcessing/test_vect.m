%file = '/archive/home/ethrane/temp/test_8359-test_835939431-60.gwf';
%file = '/archive/home/ethrane/temp/ET2/H-H1_RDS_C03_L2_ET2-835939431-30.gwf';
file = '/data/node1/ethrane/HL-SIDv1_H1L1_250mHz-816260883-26.gwf';

%channel = 'L1:LocalPSD';
%channel = 'H1L1:ImCSD';
%channel = 'H1L1:flow';
%channel = 'H1L1:fhigh';
%channel = 'H1L1:deltaF';
%channel = 'H1:AdjacentPSD';
channel = 'H1L1:w1w2bar';


gpsStart = 0;
duration = 18;
%[data,time,freq]=frgetvect(file,channel,gpsStart,duration);

gpsStart = 816260883;
segmentDuration = 52;
test = frgetvect(file, ...
			      'H1L1:w1w2bar', ...
			      gpsStart,segmentDuration)


