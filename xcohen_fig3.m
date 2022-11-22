function xcohen_fig3()

% specify frequencies
frex = linspace(3,60,50);
fwhm = linspace(.4,.1,length(frex)); % in ms
ncyc = linspace(3.2,16,length(frex));
% parameters for complex Morlet wavelets
srate = 1024;
wavtime = -2:1/srate:2;
midp = dsearchn(wavtime',0);
% outputs
empfwhm = zeros(length(frex),2);
% loop over frequencies
for fi=1:length(frex)

    % create the Gaussian using the FWHM formula (equation 3)
    gwin = exp( (-4*log(2)*wavtime.^2) ./ fwhm(fi)^2 );

    % measure the empirical fwhm
    empfwhm(fi,1) = wavtime(midp-1+dsearchn(gwin(midp:end)',.5)) - wavtime(dsearchn(gwin(1:midp)',.5));

    % create the Gaussian using the n-cycles formula (equations 1-2)
    s = ncyc(fi) / (2*pi*frex(fi));
    gwin = exp( -wavtime.^2 ./ (2*s^2) );

    % empirical FWHM
    empfwhm(fi,2) = wavtime(midp-1+dsearchn(gwin(midp:end)',.5)) - wavtime(dsearchn(gwin(1:midp)',.5));
end
figure(3), clf
plot(frex,empfwhm*1000,'o-','markersize',8,'markerfacecolor','w','linew',2)
xlabel('Wavelet frequency (Hz)'), ylabel('FWHM (ms)')
legend({'using FWHM';'Using n-cycles'})
set(gca,'xlim',[frex(1)-1 frex(end)+1])