function [z1, z2] = calc_Zs(hilbVf1,hilbVf2,hilbIf1,hilbIf2)

z1 = hilbVf1./hilbIf1;
z2 = hilbVf2./hilbIf2;

z_scale = min([max(real(z1)) - min(real(z1)), max(real(z2)) - min(real(z2))]);
eps_ = max([1e-4 * z_scale, eps]);
new_z1_real = max(real(z1),real(z2) + eps_);

factor = new_z1_real ./ real(z1);
%z1 = factor .* z1;
%z2 = factor .* z2;

filt_sz = 801;
%z1 = smoothdata(z1,"movmean",filt_sz);
%z2 = smoothdata(z2,"movmean",filt_sz);