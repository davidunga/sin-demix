function x = smooth_amplitude(x,sz)

a = smoothdata(amplitude(x),"movmean",sz);
x = set_amplitude(x,a);