function a = adjust_amp_to_new_phase(a,p,new_p)

a = a .* cos(new_p - p);
