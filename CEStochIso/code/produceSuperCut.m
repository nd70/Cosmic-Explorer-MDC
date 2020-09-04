function produceSuperCut(bt_a_5_file,bt_a0_file,bt_a3_file,outname)
%
% produceSuperCut --- function to combine bad GPS times (as identified by the delta sigma cut in makeDeltaSigmaCutsInPostProc.m) for 3 values of the spectral index alpha into a single file
%
% Input
% - bt_a_5_file = filename of delta sigma cut for a=-5
% - bt_a0_file = filename of delta sigma cut for a=0
% - bt_a3_file = filename of delta sigma cut for a=3
% - outname = filename in which to write output
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%bt_a_5=load('/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/a-5_postproc/badGPSTimes_a-5.txt'); %times of delta sigma cut for a=-5
%bt_a0=load('/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/a0_postproc/badGPSTimes_a0.txt'); %times of delta sigma cut for a=0
%bt_a3=load('/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/a3_postproc/badGPSTimes_a3.txt'); %times of delta sigma cut for a=3
%outname='/home/sgwynne.crowder/sgwb/O1_zerolag/C02/output/badTimesCut_a_5a0a3.txt'; %filename for delta sigma cut combined for a=-5,0,3
%%%

% 7/13 AAM: Added these lines
bt_a_5=load(bt_a_5_file);
bt_a0=load(bt_a0_file);
bt_a3=load(bt_a3_file);


bt_super=union(bt_a_5,union(bt_a0,bt_a3));

fid=fopen(outname,'w+');
for jj=1:length(bt_super) %write badGPSTimes to file
  fprintf(fid,'%i\n',bt_super(jj));
end
fclose(fid)
