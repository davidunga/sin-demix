function dmx = linear_demix2(v, Fs, ws, opts)

arguments
    v
    Fs
    ws
end

assert(all(ws>0));
assert(all(diff(ws) > 0), "Frequencies must be monotonically increasing");

v = v(:);
t = (0:(length(v)-1))'/Fs;

A = ones([length(t), 1+2*length(ws)]);
for k = 1:length(ws)
    A(:,2*k:2*k+1) = [sin(ws(k)*t), cos(ws(k)*t)];
end

win_dur = 2*pi/min(ws);
r = floor(win_dur * Fs / 2);  % sampling radius
sample_ixs = round(linspace(-r,r,2*r+1)); % 0-centered sampling indices

x = nan(size(A));
for i = (r+1):(size(A,1)-r)
    x(i, :) = A(i + sample_ixs, :) \ v(i + sample_ixs);
end

x(1:r,:) = repmat(x(r+1,:),[r,1]);
x(end-r+1:end,:) = repmat(x(end-r,:),[r,1]);

dc = x(:,1);
a = nan([length(ws), length(t)]);
p = nan([length(ws), length(t)]);

for k = 1:length(ws)
    i = 2*k;
    a(k,:) = sqrt(sum(x(:,i:i+1).^2,2));
    p(k,:) = unwrap(atan2(x(:,i+1), x(:,i)));
end

dmx = MixParams(Fs=Fs, ws=ws, dc=dc, a=a, p=p);
