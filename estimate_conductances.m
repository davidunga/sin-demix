function res = estimate_conductances(I,V,expr_info,cell_info,opts)

arguments
    I
    V

    expr_info                   % experiment params struct
    cell_info                   % cell params struct

    opts.demix_mode = "wt"      % demix mode
    opts.stable_demix = 0       % use stable demix
    opts.gtot_kind = 1          % 1/2. gtot computation variant. 1 produces non-negative values.
    opts.rest_dur = .2          % duration of rest window [sec] to detect. if 0, expr_info.rest_win will be used.
    opts.use_input_ws = 0       % prefer omegas from input
    opts.smooth_zs = 0          % smooth Zs
end

% ------------------------------
% Prepare:

Fs = expr_info.Fs;
ws = sort(dominant_freqs(I,Fs))*2*pi;

% validate
assert(length(ws)==2);
assert(all(ws>0));
if ~isempty(expr_info.ws)
    assert(max(abs(ws-expr_info.ws)./expr_info.ws) < .005);
end

if opts.use_input_ws
    ws = expr_info.ws;
end

period_sz = 2*pi/min(ws)*Fs;

% ------------------------------
% Demix:

Idmx = do_demix(I, Fs, ws, mode=opts.demix_mode, stable=opts.stable_demix);
Vdmx = do_demix(V, Fs, ws, mode=opts.demix_mode, stable=opts.stable_demix);

comps_hatI = Idmx.comps();
I1 = comps_hatI(2,:);
I2 = comps_hatI(3,:);

comps_hatV = Vdmx.comps();
V0 = comps_hatV(1,:);
V1 = comps_hatV(2,:);
V2 = comps_hatV(3,:);

% ------------------------------
% Impedances:

z1 = hilbert(V1)./hilbert(I1);
z2 = hilbert(V2)./hilbert(I2);

if opts.smooth_zs
    z1 = movavg(z1, period_sz);
    z2 = movavg(z2, period_sz);
end

% ------------------------------
% Electrode resistance and total conductance:

Rs = calc_Rs(cell_info.c, ws, z1, z2);

%zz1 = z1*exp(1i*.055*pi/180);
zz1=z1;
gtot = calc_gtotal(cell_info.c, ws, zz1, z2, 3);

% ------------------------------
% Clean:

Iclean = I;
for frq = ws(:)' / (2*pi)
    Iclean = bandstop(Iclean, frq + [-1, 1]*60, Fs, ImpulseResponse='fir', Steepness=.54, StopbandAttenuation=25);
end
V0 = V0 - Iclean.*real(Rs);

% ------------------------------
% Conductances:

% leak:
if opts.rest_dur == 0
    rest_ixs = get_rest_ixs_fromWindow(Fs, V0, expr_info.rest_win);
else
    rest_ixs = get_low_variance_window(V0, Fs, ws);
end
vl = mean(V0(rest_ixs));
gl = mean(gtot(rest_ixs));

% gi:
dv = diff(V0);
cdvdt = cell_info.c*[dv dv(end)]*Fs;
gi = (cdvdt + gl.*(V0 - vl) + (gtot - gl).*(V0 - cell_info.ve) - Iclean) ./ (cell_info.vi - cell_info.ve);
gi = real(gi);

% ge:
ge = gtot - gl - gi;

% ------------------------------
% Finalize:

res = struct();

res.Iin = gi.*(V0 - cell_info.vi);
res.Iex = ge.*(V0 - cell_info.ve);

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

res.Rs_est = Rs;

res.cdvdt = cdvdt;
res.Iclean = Iclean;

res.cell = cell_info;
res.expr = exprStruct(Fs=Fs,ws=ws);
return


% =================================================================
% Helper functions

function Rs = calc_Rs(c,ws,z1,z2)
q = c*(ws(2) - ws(1))*(z2 - z1);
Rs = (z1 + z2)/2 + 1i*sqrt(-q.^2 - 1i*4*q) / (2*c *(ws(2) - ws(1)));
return

function gtot = calc_gtotal(c,ws,z1,z2,kind)
% calc gtotal.
% kind = calculation variant to use. variants are mathematically equivalent, but for
% variant 1 matlab produces non-negative results.

q = c*(ws(2) - ws(1))*(z2 - z1);

%r = round(.1*length(z1));
%BIAS = .13*mean(abs(q(r:end-r)));
if kind==1
    d = q - sqrt(q.^2 + 1i*4*q);
    gtot = -2*c*(ws(2) - ws(1)) * real(1./d);
elseif kind==2
    d = q + 1i*sqrt(-q.^2 - 1i*4*q);
    gtot = -2*c*(ws(2) - ws(1)) * real(1./d);
elseif kind==3
    % approximation
    gtot = c*(ws(2) - ws(1))/sqrt(2)*abs(real((1-1i)./sqrt(q)));
else
    error("Unknown gtot kind");
end

return

% --
% Lagacy helpers:

function ixs = get_rest_ixs_fromWindow(Fs, V0, rest_win)
assert(numel(rest_win)==2 && diff(rest_win)>0);
ifm = round(rest_win(1) * Fs);
ito = round(rest_win(2) * Fs);
ixs = ifm + find(V0(ifm:ito) <= prctile(V0(ifm:ito), 15));
return

function Rs = calc_Rs_LEGACY(c,ws,z1,z2)
w1 = ws(1);
w2 = ws(2);
Rs = (1/(2*(1i*c*w1 - 1i*c*w2)))*(1i*c*w1*z1 - 1i*c*w2*z1 + 1i*c*w1*z2 - ...
    1i*c*w2*z2 + ((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 + 1i*c*w2*z2).^2 - ...
    4*(1i*c*w1 - 1i*c*w2)*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)).^0.5);
return



