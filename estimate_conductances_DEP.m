function res = estimate_conductances_DEP(I,V,rest_win,expr,cell,options)

arguments
    I
    V
    rest_win

    expr                        % experiment params struct
    cell                        % cell params struct

    options.Vrest_p = 15        % voltage percentile within rest window to classify as rest
    options.clean_r = 60        % half-size of Fpass band for cleaning
    options.demix_mode = "wt2"   % type of demix
end

% ------------------------------
% Prepare:

[w1, w2] = get_w1_w2(expr, I);

% ------------------------------
% Demix:

ws = calc_signal_properties(I, expr.Fs);

Idmx = wt_demix2c(I, expr.Fs, ws);
Vdmx = wt_demix2c(V, expr.Fs, ws);

comps_hatI = Idmx.comps();
comps_hatV = Vdmx.comps();

%[comps_hatI, Idmx] = wt_demix2(I, expr.Fs, [0, w1, w2]);
%[comps_hatI2, Idmx2] = wt_demix2b(I, expr.Fs, [0, w1, w2]);
%[comps_hatI2,~,Idmx2] = stable_demix(I, expr.Fs, w1, w2, mode=options.demix_mode, stable=true);
I1 = comps_hatI(2,:);
I2 = comps_hatI(3,:);

%[comps_hatV, Vdmx] = wt_demix2(V, expr.Fs, [0, w1, w2]);
%[comps_hatV2, Vdmx2] = wt_demix2b(V, expr.Fs, [0, w1, w2]);
%[comps_hatV2,~,Vdmx2] = stable_demix(V, expr.Fs, w1, w2, mode=options.demix_mode, stable=true);
V0 = comps_hatV(1,:);
V1 = comps_hatV(2,:);
V2 = comps_hatV(3,:);

DBG_PLT = 0;
if DBG_PLT
    t = (0:(length(V)-1))/expr.Fs;
    dmx = Vdmx;
    figure();
    tiledlayout(3,1,TileSpacing="compact",Padding="compact");
    nexttile(); hold on;
    plt x=t ax=gca dmx.a
    nexttile(); hold on;
    plt x=t ax=gca dmx.p
    nexttile(); hold on;
    plt x=t ax=gca comps_hatV(1,:)
    linkaxes(allaxes(),'x');
    xlim([.1,.9]*(t(end)-t(1)) + t(1));
end

% ------------------------------
% Impedances:

z1 = hilbert(V1)./hilbert(I1);
z2 = hilbert(V2)./hilbert(I2);

% ------------------------------
% Electrode resistance and total conductance:

Rs = calc_Rs(cell.c, w1, w2, z1, z2);
%gtot = real(1 ./ (z1 - Rs));
gtot = calc_gtotal(cell.c, w1, w2, z1, z2);

% ------------------------------
% Clean:

Iclean = I;
for frq = [w1,w2] / (2*pi)
    Fpass = frq + [-options.clean_r, options.clean_r];
    Iclean = bandstop(Iclean, Fpass, expr.Fs, ImpulseResponse='fir', Steepness=.54, StopbandAttenuation=25);
end

%V0 = V0 - Iclean .* abs(Rs);
%V0 = V0 - Iclean .* real(Rs);  % check - should be one of them

% ------------------------------
% Conductances:

% leak:
rest_ixs = get_rest_potential_ixs(expr.Fs, rest_win, V0, options.Vrest_p);
vl = mean(V0(rest_ixs));
gl = mean(gtot(rest_ixs));

% gi:
dv = diff(V0);
cdvdt = cell.c*[dv dv(end)]*expr.Fs;
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

res.Rs_est = Rs;

res.cdvdt = cdvdt;
res.Iclean = Iclean;

res.cell = cell;
res.expr = exprStruct(Fs=expr.Fs,f1=w1/(2*pi),f2=w2/(2*pi));
return


% =================================================================
% Helper functions


function ixs = get_rest_potential_ixs(Fs, rest_win, V0, Vrest_p)
ifm = round(rest_win(1) * Fs);
ito = round(rest_win(2) * Fs);
ixs = ifm + find(V0(ifm:ito) <= prctile(V0(ifm:ito), Vrest_p));
return

function [w1, w2] = get_w1_w2(expr, I)
domfrqs = dominant_freqs(I,expr.Fs,n=2);
f1 = min(domfrqs);
f2 = max(domfrqs);
if ~isnan(expr.f1)
    % omegas were given - make sure they match fft
    assert(abs(expr.f1 - f1)/expr.f1 < .01);
    assert(abs(expr.f2 - f2)/expr.f2 < .01);
    f1 = expr.f1;
    f2 = expr.f2;
end
w1 = 2 * pi * f1;
w2 = 2 * pi * f2;
return

function Rs = calc_Rs_LEGACY(c,w1,w2,z1,z2)
Rs = (1/(2*(1i*c*w1 - 1i*c*w2)))*(1i*c*w1*z1 - 1i*c*w2*z1 + 1i*c*w1*z2 - ...
    1i*c*w2*z2 + ((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 + 1i*c*w2*z2).^2 - ...
    4*(1i*c*w1 - 1i*c*w2)*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)).^0.5);
return

function Rs = calc_Rs(c,w1,w2,z1,z2)
q = c*(w2 - w1)*(z2 - z1);
Rs = (z1 + z2)/2 + 1i*sqrt(-q.^2 - 1i*4*q) / (2*c *(w2 - w1));
return

function gtot = calc_gtotal(c,w1,w2,z1,z2)
q = c*(w2 - w1)*(z2 - z1);
gtot = -2*c*(w2 - w1) * real(1./(q + 1i*sqrt(-q.^2 - 1i*4*q)));
return

function gtot = calc_gtotal_variant1(c,w1,w2,z1,z2)
q = c*(w2 - w1)*(z2 - z1);
gtot=-2*c*(w2 - w1) * real(1./(q - sqrt(q.^2 + 1i*4*q)));
return
