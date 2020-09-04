function [] = bad()

% E. Thrane
% Run before running sumXFisher_S5.  Identify all segments failing the sigma
% cut in the bad_segs file.  If there is an error loading the naivesigmas file 
% for a job, the entire job is listed in the bad_jobs file.  Example of a bad
% job: job 6 runs w/o errors on stochastic, but there are no segments in the
% output file.

bad_segs=fopen('/archive/home/ethrane/SID/bad_segs.txt','w+');
bad_jobs=fopen('/archive/home/ethrane/SID/bad_jobs.txt','w+');
for job=1:18837
  gpsfile = ['~/S5H1L1_sph52_SID_L20/S5H1L1_sph52_SID_gps_naivesigmas.job' ... 
             num2str(job) '.trial1.dat'];
  try
    gpsData = load(gpsfile);
    start = gpsData(:,1);
    naive = gpsData(:,2);
    theor = gpsData(:,3);

    for ii=1:length(start)
      if theor(ii)/naive(ii)<0.8 || theor(ii)/naive(ii)>1.2  % sigma cut
        fprintf(bad_segs,'%i %i %i\n',job,ii,start(ii));
      end
    end
  catch
    fprintf(bad_jobs,'%i\n',job);
  end
end  
fclose(bad_segs);
fclose(bad_jobs);
