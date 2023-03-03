function res = estimate_conductances(I,V,rest_win,expr,cell,options)

arguments
    I
    V
    rest_win

    expr                        % experiment params struct
    cell                        % cell params struct

    options.Vrest_p = 15        % voltage percentile within rest window to classify as rest
    options.clean_r = 60        % half-size of Fpass band for cleaning
    options.demix_mode = "wt"   % type of demix
end

% ------------------------------
% Prepare:

[w1, w2] = get_w1_w2(expr, I);
Fs = expr.Fs;

% ------------------------------
% Demix:

[comps_hatV,~] = stable_demix(V, Fs, w1, w2, mode=options.demix_mode);
V0 = comps_hatV(1,:);
V1 = comps_hatV(2,:);
V2 = comps_hatV(3,:);

[comps_hatI,~] = stable_demix(I, Fs, w1, w2, mode=options.demix_mode);
I1 = comps_hatI(2,:);
I2 = comps_hatI(3,:);

% ------------------------------
% Impedances:

hV1 = hilbert(V1);
hI1 = hilbert(I1);
hV2 = hilbert(V2);
hI2 = hilbert(I2);
z1 = hV1./hI1;
z2 = hV2./hI2;

% ------------------------------
% Electrode resistance and total conductance:

Rs_est = calc_Rs(cell.c, w1, w2, z1, z2);
gtot = 1i*(1i + cell.c*Rs_est*w1 - cell.c*w1*z1) ./ (Rs_est - z1);
gtot = real(gtot);

% ------------------------------
% Clean:

Iclean = I;
for frq = [w1,w2] / (2*pi)
    Fpass = frq + [-1, 1] * options.clean_r;
    Iclean = bandstop(Iclean, Fpass, Fs, ImpulseResponse='fir', Steepness=.54, StopbandAttenuation=25);
end

V0 = V0 - Iclean .* abs(Rs_est);
V0 = V0 - Iclean .* real(Rs_est);

% ------------------------------
% Conductances:

% leak:
rest_ixs = get_rest_potential_ixs(Fs, rest_win, V0, options.Vrest_p);
vl = mean(V0(rest_ixs));
gl = mean(gtot(rest_ixs));

% gi:
dv = diff(V0);
cdvdt = cell.c*[dv dv(end)]*Fs;
gi = (cdvdt + gl.*(V0 - vl) + (gtot - gl).*(V0 - cell.ve) - Iclean) ./ (cell.vi - cell.ve);
gi = real(gi);

% ge:
ge = gtot - gl - gi;

% ------------------------------
% Finalize:

res = struct();

res.Iin = gi.*(V0 - cell.vi);
res.Iex = ge.*(V0 - cell.ve);

res.ge = ge;
res.gi = gi;
res.gl = gl;
res.gtot = gtot;

res.z1 = z1;
res.z2 = z2;

res.V0 = V0;
res.V1 = V1;
res.V2 = V2;

res.I1 = I1;
res.I2 = I2;

res.Rs_est = Rs_est;

res.cdvdt = cdvdt;
res.Iclean = Iclean;

res.cell = cell;
res.expr = exprStruct(Fs=Fs,f1=w1/(2*pi),f2=w2/(2*pi));
return


% =================================================================
% Helper functions


function ixs = get_rest_potential_ixs(Fs, rest_win, V0, Vrest_p)
ifm = round(rest_win(1) * Fs);
ito = round(rest_win(2) * Fs);
ixs = ifm + find(V0(ifm:ito) <= prctile(V0(ifm:ito), Vrest_p));
return

function [w1, w2] = get_w1_w2(expr, I)
[FT, FT_frqs] = fourier(I, expr.Fs);
[~, dominant_frqs] = findpeaks(abs(FT), FT_frqs, MinPeakDistance=5, NPeaks=2);
w1 = min(dominant_frqs) * 2 * pi;
w2 = max(dominant_frqs) * 2 * pi;
if expr.w1 ~= 0
    % omegas were given - make sure they match fft
    assert(abs(expr.w1 - w1)/expr.w1 < .01);
    assert(abs(expr.w2 - w2)/expr.w2 < .01);
    w1 = expr.w1;
    w2 = expr.w2;
end
return

function Rs = calc_Rs(c,w1,w2,z1,z2)
Rs = (1/(2*(1i*c*w1 - 1i*c*w2)))*(1i*c*w1*z1 - 1i*c*w2*z1 + 1i*c*w1*z2 - ...
    1i*c*w2*z2 + ((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 + 1i*c*w2*z2).^2 - ...
    4*(1i*c*w1 - 1i*c*w2)*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)).^0.5);
return
