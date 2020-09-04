function [cal] = getTrivialCalibration()
%
%  produce a dummy callibration
%
%  Routine written by Stefan Ballmer
%
%  $Id: getTrivialCalibration.m,v 1.1 2007-12-11 16:54:01 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cal.t=[0;999999999];
cal.f=(0:0.125:7000)';
cal.R0=ones(size(cal.f));
cal.C0=cal.R0;
cal.alpha=ones(size(cal.t));
cal.gamma=ones(size(cal.t));

return
