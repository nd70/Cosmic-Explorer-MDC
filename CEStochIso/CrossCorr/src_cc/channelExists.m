function result=channelExists(channelName, frameFile)
%
% channelExists(channelName, frameFile) - returns true if the channel
% exists in the specified frame file, otherwise false
%
% - channelName is a string containing the channel name
% - frameFile is the path to a frame file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  result = true;

  % Before we look for the channel we need to make sure
  % the frame file exists, since frgetvect only gives a
  % warning
  checkFileExists('Frame file', frameFile);

  % Store the last warning so we can restore it later. This
  % is needed in case users rely on lastwarn in other parts
  % of their code
  [lastMsg, lastWarn] = lastwarn;

  % Set the last warning to the empty string
  warning('');

  % This will read in 1 second of data from the start of the frame and
  % give a warning if the channel is not found.
  % We really need a MEX file which can check for the existence
  % of a channel more efficiently
  data = frgetvect(frameFile, channelName);
  [frMsg, frWarn] = lastwarn;

  % Silently restore the original last warning
  lastwarn(lastMsg, lastWarn);

  if (strcmp(frWarn, 'frgetvect:channelNotFound'))
    result = false;
  end;

return;
