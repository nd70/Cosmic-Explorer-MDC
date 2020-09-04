%file = '/archive/home/ethrane/temp/ET2/H-H1_RDS_C03_L2_ET2-835939491-30.gwf';
file = '/data/node1/ethrane/HL-SIDv1_H1L1_250mHz-816260883-26.gwf';

%channel = 'ParametersChannelName';
channel = 'H1L1:Params';

params=frgetparams(file,channel);
