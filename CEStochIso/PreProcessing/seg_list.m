% Create a list of all segment start times based on naive sigma files.
% Calculate ratio for sigma cut while were at it.  -E. Thrane
%
clear all;
NJobs = 18837;    % total number of jobs to combine
segs.gps = [];
segs.job = [];
segs.NFail = 0;
for ii=1:NJobs
  try
    naive_dat = load(['/archive/home/ethrane/preproc/naive/' ...
     'S5H1L1_sph52_noIM_mask5_L12_naivesigmas.job' num2str(ii) '.trial1.dat']);
    segs.gps = [segs.gps ; naive_dat(:,1)];
    segs.job = [segs.job ; ii*ones(size(naive_dat(:,1)))];
  catch
    segs.NFail = segs.NFail + 1;
    fprintf('Error with job %i.\n',ii);
  end
end

save('segs.mat','segs');
