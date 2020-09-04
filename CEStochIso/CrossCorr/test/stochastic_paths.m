% Eric Thrane
%
% This is an example of a startup.m file which you can run automatically 
% every time you launch matlab to set your matlab paths.  These commands will
% set up the paths needed for stoch code.  You must define the stoch_install
% variable to point to your svn copy of stamp.
% pwd uses the home directory.  if you want to hard-code a path, that's fine
% too.  Just un-comment out the subsequent line and modify for your own home
% directory.

%stoch_install = pwd;
%stoch_install = '/home/ethrane/matapps/';
matapps_dir = '/home/charlton/Analysis/matapps';
stoch_install = '/home/charlton/Analysis/matapps/packages/stochastic/trunk';
stoch_tag = '/home/charlton/Analysis/matapps/packages/stochastic/tags/stochastic_v3';

%------------------------------------------------------------------------------

% make sure that the stamp install exists
if ~exist([stoch_install '/stamp2/src'])
  fprintf('stoch_install = %s is not the correct directory.\n', stoch_install);
  fprintf('Please modify stochastic_paths.m\n');
end

% define STAMP paths
% include everything in the src/ directory including subdirectories
fprintf('loading stochastic packages...');
addpath(genpath([stoch_install '/stamp2/src']));
addpath(genpath([stoch_install '/stamp2/input']));
addpath([ stoch_install '/Utilities' ]);

addpath([ matapps_dir '/admin/utilities/Channel' ]);
addpath([ matapps_dir '/admin/utilities/detgeom/matlab' ]);
addpath([ matapps_dir '/admin/utilities/misc' ]);
addpath([ matapps_dir '/admin/utilities/not4toolbox/FTSeries' ]);

% Add path to files modified from trunk last so they go at the start of the path
addpath(genpath([ stoch_tag '/CrossCorr/src_cc' ]));
addpath(genpath([ stoch_tag '/PostProcessing' ]));

fprintf('done.\n');
