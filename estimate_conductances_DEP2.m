function res = estimate_conductances(I,V,expr,cell,opts)

arguments
    I
    V

    expr                        % experiment params struct
    cell                        % cell params struct

    opts.demix_mode = "wi"        % demix mode
    opts.stable_dmx = true      % use stable demix
    opts.fast_sp = false        % FAST signal properties estimation
    opts.clean_r = 60           % half-size of Fpass band for cleaning
    opts.gtot_kind = 2          % 1/2 - function to calculate gtot. kind=2 produces non-negative values.
    opts.rest_dur = .2          % duration of rest window [sec] to detect
end

% ------------------------------
% Prepare:

% calc omegas, margin of influence, etc.
if ~opts.fast_sp
    spI = calc_signal_properties_FAST(I, expr.Fs);
    spV = spI;
else
    spI = calc_signal_properties(I, expr.Fs);
    spV = calc_signal_properties(V, expr.Fs);
end

% validate omegas
assert(length(spV.ws)==2); % we expect two omegas
assert(max(abs(spV.ws-spI.ws)./spI.ws) < .005);  % omegas computed from I and V should be equal (upto tolerance)
if ~isempty(expr.ws)
    % if omegas were given - check they are equal to computed, upto tolerance
    assert(max(abs(spV.ws-expr.ws)./expr.ws) < .005);
end

w1 = (spV.ws(1)+spI.ws(1))/2;
w2 = (spV.ws(2)+spI.ws(2))/2;

% ------------------------------
% Demix:

dmx_paramsI = {};
dmx_paramsV = {};
if string(opts.demix_mode) == "wi" && ~opts.fast_sp
    % avoid recomputing cwt..
    dmx_paramsI = {'cwtRes', spI.cwtRes};
    dmx_paramsV = {'cwtRes', spV.cwtRes};
end

Idmx = stable_demix2(I, expr.Fs, [w1,w2], demix_mode=opts.demix_mode, demix_params=dmx_paramsI, stable=opts.stable_dmx);
Vdmx = stable_demix2(V, expr.Fs, [w1,w2], demix_mode=opts.demix_mode, demix_params=dmx_paramsV, stable=opts.stable_dmx);

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

% ------------------------------
% Electrode resistance and total conductance:

Rs = calc_Rs(cell.c, w1, w2, z1, z2);
gtot = calc_gtotal(cell.c, w1, w2, z1, z2, opts.gtot_kind);

% ------------------------------
% Clean:

Iclean = I;
for frq = [w1,w2] / (2*pi)
    Fpass = frq + [-opts.clean_r, opts.clean_r];
    Iclean = bandstop(Iclean, Fpass, expr.Fs, ImpulseResponse='fir', Steepness=.54, StopbandAttenuation=25);
end
V0 = V0 - Iclean.*real(Rs);

% ------------------------------
% Conductances:

% leak:
rest_ixs = calc_rest_ixs(expr.Fs, V0, opts.rest_dur, spV.margin);
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
res.expr = exprStruct(Fs=expr.Fs,ws=[w1,w2]);
return


% =================================================================
% Helper functions


function ixs = calc_rest_ixs(Fs, V0, rest_dur, margin)
xx = movstd(V0,round(rest_dur*Fs));
[~,mni] = min(xx(margin:end-margin));
r = round(rest_dur*Fs/2);
ixs = (mni+margin-r):(mni+margin+r);
return


function ixs = get_rest_ixs_LEGACY(Fs, V0, rest_win)
ifm = round(rest_win(1) * Fs);
ito = round(rest_win(2) * Fs);
ixs = ifm + find(V0(ifm:ito) <= prctile(V0(ifm:ito), 15));
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

function gtot = calc_gtotal(c,w1,w2,z1,z2,mode)
assert(any(mode==[1,2]));
q = c*(w2 - w1)*(z2 - z1);
if mode==1
    gtot = -2*c*(w2 - w1) * real(1./(q + 1i*sqrt(-q.^2 - 1i*4*q)));
elseif mode==2
    gtot=-2*c*(w2 - w1) * real(1./(q - sqrt(q.^2 + 1i*4*q)));
end
return

