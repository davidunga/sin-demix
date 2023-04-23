% dev01

f1=100;
f2=200;
dur = 2;
Fs = 5000;

[a1,t] = make_test_signal(Fs, dur, 0);
[a2,~] = make_test_signal(Fs, dur, .2);

amp_max = .02;
amp_min = .8 * amp_max;
a1 = amp_min + (amp_max - amp_min) * a1;
a2 = amp_min + (amp_max - amp_min) * a2;

dc_min = -.07;
dc_max = -.03;
[dc,~] = make_test_signal(Fs, dur, .1);
dc = dc_min + (dc_max - dc_min) * dc;


p1 = 0;
p2 = 0;
s1 = a1 .* sin(2*pi*f1*t + p1);
s2 = a2 .* sin(2*pi*f2*t + p2);
v = dc + s1 + s2;
