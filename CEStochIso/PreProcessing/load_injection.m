function [h1 h2] = load_injection(params, stamp_inj, fs)
% function [h1 h2] = load_injection(params, stamp_inj, fs)
% e.g., [h1 h2] = load_injection(params, stamp_inj, 4096);
% params = parameter struct
% stamp_inj = stamp parameter struct
% fs = sampling frequency
% If the injection is of type 'PSD', this code loads in the f, PSD(f) data.
% If the injection is of type 'time', this code loads in t, hp(t), hx(t) data.
% by E. Thrane

  try
    stamp_inj.file;
    stamp_inj.type;
    stamp_inj.ra;
    stamp_inj.decl;
    stamp_inj.start;
  catch
    error('Error: missing variables in stamp_inj struct.\n');
  end

  % get detector information
  site1 = getsitefromletter(params.ifo1(1));
  site2 = getsitefromletter(params.ifo2(1));
  det1 = getdetector(site1);
  det2 = getdetector(site2);

  % fill pp struct
  pp.flow = params.flow;
  pp.fhigh = params.fhigh;
  pp.deltaF = params.deltaF;

  switch stamp_inj.type
    case 'PSD';
      % inj_file is of the form: freq, PSD(freq)
      Hf = load(stamp_inj.file);
    case 'time_series';
      % inj_file is of the form: time, hplus, hcross
      Hf = load(stamp_inj.file);
      % make sure signal is divisible by sampling time
      if mod(Hf(end,1), 1/fs)~=0
        Hf = [Hf ; Hf(end,1)-mod(Hf(end,1), 1/fs)+1/fs 0 0];
      end
  end

  % Set defaults for getPointSourceData_IM
  amplitude = 1;       % do not scale the injection
  intLog = 1;          % interpolate data is always on
  MakeIncoherent = 0;  % coherent source
 
  % get the whole entire series for the injection duration
  [h1, h2] = getPointSourceData_IM(stamp_inj.type, ...
    stamp_inj.start, stamp_inj.dur, pp, fs, Hf, det1, det2, ...
    stamp_inj.ra, stamp_inj.decl, amplitude, MakeIncoherent, intLog, params);

return
