function x = dc_smooth(x,win,opt)
% smooth dc component of signal

arguments
    x
    win
    opt.method = "gaussian"
end

[a,u,l,s,ui,li] = amplitude(x,1);
x = x - s + smoothdata(s,opt.method,win);