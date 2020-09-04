function [startGPS, endGPS, params] = check_times(params)

  % Check if start/end GPS times are defined. ---------------------------------
  % If not, uses the start/end times from the jobfile.
  if (isfield(params,'startGPS'))
    if (params.startGPS >= params.pass.GPStimes(1))
      startGPS = params.startGPS;
    else
      if (params.startGPS < params.pass.GPStimes(1))
        warning('The startGPS time specified in the paramfile is earlier than that provided by the jobfile... using the startGPS time specified by the jobfile.');
      end
      startGPS = params.pass.GPStimes(1);
    end
      else
	startGPS = params.pass.GPStimes(1);
  end

  which_beg_pos=min(abs(params.pass.GPStimes-startGPS));
  beg_pos = min(find(abs(params.pass.GPStimes-startGPS)==which_beg_pos));

  if (isfield(params,'endGPS'))
    if (params.endGPS <= params.pass.GPStimes(end))
      endGPS = params.endGPS;
    else
      if (params.endGPS > params.pass.GPStimes(end))
      warning('The endGPS time specified in the paramfile is later than that provided by the jobfile... using the endGPS time specified by the jobfile.');
      end
      endGPS = params.pass.GPStimes(end);
    end
  else
    endGPS = params.pass.GPStimes(end);
  end

  which_end_pos=min(abs(params.pass.GPStimes-endGPS));
  end_pos = max(find(abs(params.pass.GPStimes-endGPS)==which_end_pos));

  params.pass.which_segs=beg_pos:1:end_pos;
  params.pass.GPStimes=params.pass.GPStimes(params.pass.which_segs);
