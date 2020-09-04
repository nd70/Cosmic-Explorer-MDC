% Script to test post-processing of data from example_server_all output

doRenormalize = false;
modifyFilter = false;
displayResults = false;
applyBadGPSTimes = false;
badGPSTimesFile = 'badGPSTimes.dat';
plotResults = false;
printFormats = {};

if (exist('example_resumeCRfMJ.mat'))
  delete('example_resumeCRfMJ.mat');
end;

if (exist(badGPSTimesFile))
  delete(badGPSTimesFile);
end;

if (exist('test_postProcess.out'))
  delete('test_postProcess.out');
end;
diary('test_postProcess.out');
postProcessScriptFull('example_server_all.txt', 'example_jobs.txt', '.', 0.2, 1, ...
                      doRenormalize, modifyFilter, displayResults, applyBadGPSTimes, ...
                      badGPSTimesFile, printFormats);
diary off;
