% clear all;
% designed to work well in S5 job 3
%
%

fs = 16384;
f0 = 90;
dur = 300;
%A = 7e-22;
A = 1.5e-22;
rampdur = 2;

t=[1/fs : 1/fs: dur];
ramp = ones(size(t));
ramp(1:rampdur*fs+1) = ramp(1:rampdur*fs+1) .* [0:1/(fs*rampdur):1];
ramp(end-rampdur*fs:end) = ramp(end-rampdur*fs:end) .* [1:-1/(fs*rampdur):0];

hp = ramp.*A.*sin(2*pi*f0*t);
hx = ramp.*A.*sin(2*pi*f0*t+pi/2);

out = [t ; hp ; hx]';
save wave90Hz.dat out -ASCII;


