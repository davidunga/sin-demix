function [z1, z2, hI1, hI2, hV1, hV2] = calc_Zs(I1,I2,V1,V2)

hV1 = hilbert(V1);
hI1 = hilbert(I1);
hV2 = hilbert(V2);
hI2 = hilbert(I2);
z1 = hV1./hI1;
z2 = hV2./hI2;

z_scale = min([max(real(z1)) - min(real(z1)), max(real(z2)) - min(real(z2))]);
eps_ = max([1e-4 * z_scale, eps]);
new_z1_real = max(real(z1),real(z2) + eps_);

factor = new_z1_real ./ real(z1);
%z1 = factor .* z1;
%z2 = factor .* z2;

filt_sz = 501;
z1 = smoothdata(z1,"movmean",filt_sz);
z2 = smoothdata(z2,"movmean",filt_sz);