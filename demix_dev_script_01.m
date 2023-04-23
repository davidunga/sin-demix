% demix_dev_script_01

Fs = 1000;
dur = 3;
N = round(dur * Fs);
t = (0:(N-1))/Fs;

ws = [0, 100, 200];

dc = piecewise_const(N, [0, 2, .5], [0, .25, .75]);
amp1 = piecewise_const(N, [5, 4.5], [0, .4]);
amp2 = piecewise_const(N, [6.5], [1]);
ph1 = -.1 * ones(size(t));
ph2 = pi/8 * ones(size(t));
s1 = amp1 .* sin(ws(2)*t + ph1);
s2 = amp2 .* sin(ws(3)*t + ph2);

v = dc + s1 + s2;

% -----

[comps, opt] = wt_demix2(v, Fs, ws);

%[comps, opt] = linear_demix_fixedphase(v, Fs, ws(ws>0), opt.ph);

figure();
hold on;
plot(t, dc, 'r');
plot(t, amp1, 'r');
plot(t, amp2, 'm');
plot(t, comps(1,:), 'b');
plot(t, opt.a(1,:), 'b');
plot(t, opt.a(2,:), 'c');
title('Amplitude');

figure();
hold on;
plot(t, ph1, 'r');
plot(t, ph2, 'm');
plot(t, opt.p(1,:), 'b');
plot(t, opt.p(2,:), 'c');
plot(t, opt.ph(1)*ones(size(t)), ':b');
plot(t, opt.ph(2)*ones(size(t)), ':c');
title('Phase');

%disp(opt.ws_refined);

%plt x=t comps(1,:) comps(2,:) dc s1

% [comps, opt] = wt_demix2(v, Fs, ws, interp_f=true, interp_t=true);
% plt x=t opt.a(1,:) opt.a(2,:) opt.a(3,:) dc
% title("interp_f=true, interp_t=true");
% 
% [comps, opt2] = wt_demix2(v, Fs, ws, interp_f=false, interp_t=false);
% plt x=t opt2.a(1,:) opt2.a(2,:) opt2.a(3,:) dc
% title("interp_f=false, interp_t=false");
% 
% 
% plt x=t opt.a(1,:) opt.a(2,:) opt.a(3,:) opt2.a(1,:) opt2.a(2,:) opt2.a(3,:)
