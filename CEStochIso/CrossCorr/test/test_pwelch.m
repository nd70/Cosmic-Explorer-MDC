function test_pwelch()
% Script to verify that pwelch() will reproduce the results of psd()

  warning('off', 'signal:psd:PSDisObsolete');

  global TOL detrendFlag srate dt;
  TOL = 1.0E-12;
  detrendFlag = 'none';
  srate = 16384;
  dt = 1/srate;

  pass = true;

  pass = pass & test_real(2, 'boxcar');
  pass = pass & test_real(256, 'boxcar');
  pass = pass & test_real(1024, 'boxcar');

  % Require minimum window size of 4 for a sensible Hann-windowed psd
  pass = pass & test_real(4, 'hann');
  pass = pass & test_real(256, 'hann');
  pass = pass & test_real(1024, 'hann');

  pass = pass & test_complex(2, 'boxcar');
  pass = pass & test_complex(128, 'boxcar');
  pass = pass & test_complex(512, 'boxcar');

  % Require minimum window size of 4 for a sensible Hann-windowed psd
  pass = pass & test_complex(4, 'hann');
  pass = pass & test_complex(128, 'hann');
  pass = pass & test_complex(512, 'hann');

  pass = pass & test_sides(512, 'hann');

  if (pass)
    fprintf('PASS: all\n');
  else
    fprintf('FAIL\n');
  end;

return;

function pass=test_real(fftLength, wStr)

  global TOL detrendFlag srate dt;

  dataLength = 12*fftLength;
  overlap = fftLength/2;

  fprintf('** Testing real data length %d, %s window.\n', fftLength, wStr);

  pass = true;
  w = str2func(wStr);
  window = w(fftLength);
  data = randn(dataLength, 1);

  [psdOut, psdFreq] = psd(data, fftLength, srate, window, overlap, detrendFlag);
  
  [pwOut, pwFreq] = pwelch(data, window, overlap, fftLength, srate);

  % PWELCH output is 2*dt*PSD output for n = 2..N/2 and dt*PSD for n = 1, N/2+1
  pwOut(2:end-1) = 1/(2*dt)*pwOut(2:end-1);
  pwOut(1) = (1/dt)*pwOut(1);
  pwOut(end) = (1/dt)*pwOut(end);

  m = max(abs(psdOut - pwOut));
  if (m < TOL)
    fprintf('PASS')
  else
    pass = false;
    fprintf('FAIL');
  end;
  fprintf(' PSD values with max. difference %f\n', m);

  m = max(abs(psdFreq - pwFreq));
  if (m < TOL)
    fprintf('PASS')
  else
    pass = false;
    fprintf('FAIL');
  end;
  fprintf(' frequency with max. difference %f\n', m);

return;

function pass=test_complex(fftLength, wStr)

  global TOL detrendFlag srate dt;

  dataLength = 12*fftLength;
  overlap = fftLength/2;

  fprintf('** Testing complex data length %d, %s window.\n', fftLength, wStr);

  pass = true;
  w = str2func(wStr);
  window = w(fftLength);
  data = randn(dataLength, 1) + i*randn(dataLength, 1);

  [psdOut, psdFreq] = psd(data, fftLength, srate, window, overlap, detrendFlag);
  
  [pwOut, pwFreq] = pwelch(data, window, overlap, fftLength, srate);

  % PWELCH output is dt*PSD output for all n
  m = max(abs(dt*psdOut - pwOut));
  if (m < TOL)
    fprintf('PASS')
  else
    pass = false;
    fprintf('FAIL');
  end;
  fprintf(' PSD values with max. difference %f\n', m);

  m = max(abs(psdFreq - pwFreq));
  if (m < TOL)
    fprintf('PASS')
  else
    pass = false;
    fprintf('FAIL');
  end;
  fprintf(' frequency with max. difference %f\n', m);

return;

function pass=test_sides(fftLength, wStr)

  global TOL detrendFlag srate dt;

  dataLength = 12*fftLength;
  overlap = fftLength/2;

  fprintf('** Testing one-sided vs. two-sided real data length %d, %s window.\n', fftLength, wStr);

  pass = true;
  w = str2func(wStr);
  window = w(fftLength);
  data = randn(dataLength, 1);

  [pwOut1, pwFreq] = pwelch(data, window, overlap, fftLength, srate);
  [pwOut2, pwFreq] = pwelch(data, window, overlap, fftLength, srate, 'twosided');

  % One-sided psd is twice the two-sided psd except at DC (1) and Nyquist (N/2+1)
  pwOut2(2:fftLength/2) = 2*pwOut2(2:fftLength/2);

  m = max(abs(pwOut1 - pwOut2(1:length(pwOut1))));
  if (m < TOL)
    fprintf('PASS')
  else
    pass = false;
    fprintf('FAIL');
  end;
  fprintf(' PSD values with max. difference %f\n', m);

return;
