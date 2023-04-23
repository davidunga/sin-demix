function s = mixParamsStruct(params)

arguments
    params.Fs
    params.ws
    params.dc
    params.a
    params.p = []
    params.ph = []
end

Fs = params.Fs;
ws = params.ws;
a = params.a;
p = params.p;
dc = params.dc;
ph = params.ph;

if iscell(a)
    a = cell2mat(a(:));
end
if isempty(p)
    assert(~isempty(ph))
    p = repmat(ph(:),[1,size(a,2)]);
end
if iscell(p)
    p = cell2mat(p(:));
end

dc = dc(:)';
ws = ws(:)';

assert(size(a,1)==size(p,1));
assert(size(a,2)==size(p,2));
assert(length(dc)==size(a,2));
assert(length(ws)==size(a,1));
assert(isempty(ph) || length(ph)==length(ws));

s = struct();
s.Fs = Fs;
s.ws = ws;
s.dc = dc;
s.a = a;
s.p = p;
s.ph = ph;
s.t = (0:(length(dc)-1))/Fs;
