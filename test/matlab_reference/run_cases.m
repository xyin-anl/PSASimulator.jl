% Evaluate rows of cases.csv with the original MATLAB code (PEESEgroup/PSA).
% Run from a checkout of https://github.com/PEESEgroup/PSA:
%   octave-cli run_cases.m <cases.csv> <first_row> <last_row> <out.csv>
args = argv(); cases = args{1}; r0 = str2double(args{2}); r1 = str2double(args{3}); out = args{4};
addpath('CycleSteps'); load('Params');
C = csvread(cases); N = 10;
fid = fopen(out, 'w');
for r = r0:r1
  s = C(r,1); idx = C(r,2); x = C(r,4:9);
  if C(r,3) == 0, typ = 'ProcessEvaluation'; else, typ = 'EconomicEvaluation'; end
  material = {SimParam(idx,:), IsothermPar(idx,:)};
  % Same mapping as PSACycleSimulation.m, without its try/catch
  pv = [1.0, x(1), x(1)*x(4)/8.314/313.15, x(2), x(3), x(5), 1e4, x(6)];
  tic;
  [o, c] = PSACycle(pv, material, [], typ, N);
  fprintf(fid, '%d,%d,%.15g,%.15g,%.15g,%.15g,%.15g,%.2f\n', s, idx, o(1), o(2), c(1), c(2), c(3), toc);
  fflush(fid);
end
fclose(fid);
