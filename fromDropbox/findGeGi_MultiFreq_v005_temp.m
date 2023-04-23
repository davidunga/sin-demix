function [ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex] = findGeGi_MultiFreq_v005_temp(V,I,sf,c,reves,searchtime,plotit,FILTP,cValue,BoostCe,cableBoost,FreqArray,FiltType,hybridCe,removeHighLow,compensLowRin);
%simple way to find conductances.
%based on injection of TWO frequencies (for example 400 and 900 Hz).
% the lowest freq is used for calculation of Ge and Gi and the higher one
% to calculate the Rs as we assume shortcut of the cell at this frequency.
%cValue of <0  will caluclate c automatically
% dec 6 2020 adding solution for 4 different frequencies as a way to find also
%the capcitance of the electrode (i.e., re,ce, c_cell, g_cell).
% cableBoost 1 for assuming a cable structure, can help a lot
%FreqArray May 2021 in case of more than two frequencies, select the
%combination for calculation. e.g., 144, 190, 230 270 Hz, [1 3] will use
%144 and 230 Hz. This allows averaging of the combination of Ge Gi for when
%measurfements are noisy.
% hybridCe added 14/2/2022 to calculate c with the assumption of pippete c
% but with original z1 and Z2
% last 24_02_2022 4:33 pm
%version 005 with removal of artifact by bandpass the voltage
%at higher and lower frequencies than ff and ff2 and
%subtraction of the filtered voltage from V.
%last Aug 7 2022
global evf1 evf2 VF1 VF2 IF1 IF2 IF eif1 eif2 hilbVf1 hilbIf1 hilbVf2 hilbIf2 ff ff2 FFTV VC
global DFilter


VVo = V;

Xfound = [0 0];

global GTC;
sff = sf
global cmm
ge = 0;
gi = 0;
Rs_meas= 0;
if ~exist('plotit');

    plotit = 0;
end

if ~exist('compensLowRin');
    compensLowRin = 0;
end




cmm = cValue;

modec = 1;

if exist('cValue'); % if cValue>0 c is not automaticly calculated and taken from this value
    if cValue > 0
        modec = 0;
        c = cValue;
    end
end

if exist('FiltType'); %base on 'bandpass' (==1) or 'filtfilt' (==0)
    Filtype = FiltType; % 0 or 1
else
    Filtype = 1; % default on bandpass
end

if isempty(c)% calculate c from the sinosiodal data based on some assumptions that:
    %Z = R+g/(g^2+w^2*g^2)-j*w*c/(g^2+w^2*g^2); R electrode resist
    % if w*c>>g than
    %Z ~= R-j*(1/w*c);
    modec = 1;

end

hyBridCe = 0;
if exist('hybridCe');
    hyBridCe = hybridCe;
end

%searching for the two frequencies:
MC = 1000000;
autofreq= 1;

dt = 1/sf;
sf11 = sf;

df = 1./(dt*length(V));
fv = abs(fft(I-mean(I)));
FFTV = fv;
maxfv = max(fv(round(100/df):end));

fvshort = fv(round(100/df):end/2);
factmaxstd = maxfv/std(fvshort);
global MPH;
MPH = maxfv*0.55;
%MPH = std(fv)*5;

[pl,lc] = findpeaks(fv(round(100/df):round(end/2.2)),'MinPeakHeight',MPH);
[pl,lc] = findpeaks(fv(round(100/df):end-round(100/df)),'MinPeakHeight',MPH);
%[pl,lc] = findpeaks(fv(round(100/df):(1300/dt)),'MinPeakHeight',MPH);
global AAAA
global LCC
LCC = lc;
AAAA = fvshort;
%remove near peaks:
templc = [];;
wfr = round(df/10); %window points


if ~exist('FreqArray'); %May 2021 allows to select the frequenies.

    ff = lc(1)*df+99;
    ff2 = lc(2)*df+99;
else

    ff = lc(FreqArray(1))*df+99;
    ff2 = lc(FreqArray(2))*df+99;
end


'done lc'


%%%
%%% new 2022 procssing of VC and subtract it before other things
VCfirst = 0;
if VCfirst % processing VC
    VA = V;
    filyert = 'fir';
    DFF = 65;
    A= 20;
    VC = V;
    Iclean = I;
    LCC = lc;
    if 1 % essential do not change

        for MF = 1:length(lc)/2;
            ff3toend = lc(MF)*df+99;
            STF3 = 0.62;
            VC = bandstop(VC,[ff3toend-DFF ff3toend+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',A);
            Iclean = bandstop(Iclean,[ff3toend-DFF ff3toend+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',A);
            global VC1
            VC1 = VC;
        end
    end
    %return
    %lowpass the VC and Iclean below the lower frequency of the two
    %frequencues Jan 2022;

    if 0
        if Filtype == 1 %when using 'filt' to filter V
            VC =  lowpass(VC,ff-10,1/dt,'ImpulseResponse',filyert,'Steepness',0.997,'StopbandAttenuation',A);
            Iclean =  lowpass(Iclean,ff-10,1/dt,'ImpulseResponse',filyert,'Steepness',0.9997,'StopbandAttenuation',A) ;
            VC =  lowpass(V,ff-30,1/dt,'ImpulseResponse','iir','Steepness',0.999,'StopbandAttenuation',120);
            Iclean =  lowpass(I,ff-30,1/dt,'ImpulseResponse','iir','Steepness',0.999,'StopbandAttenuation',120) ;
            global VC2
            VC2 = VC;

        else % using the 'filt' option
            bpFilt_FFVclean=  designfilt('lowpassiir', 'FilterOrder', 10, ...
                'PassbandFrequency', ff-10, 'PassbandRipple', 0.2,...
                'SampleRate', 1/dt);
            VC = filtfilt(bpFilt_FFVclean,VC);
            Iclean = filtfilt(bpFilt_FFVclean,Iclean);
            'lastone....'
        end

    end
    V = V-VC; %Sep 2022 to remove VC before computing filters.
end

%% end VC processing.

%% finding envelopes:
DFF = 12.0; %for bonn 22
dst = 0.11; %%for bonn 0.16
NC = 3; %2 for bonn 2
STF1 = 0.57; %larger value smeers the time course!!! for bonn 0.57;
STF2 = 0.57;
if nargin>7
    if exist('FILTP');

        DFF = FILTP(1); %for bonn 22
        dst = FILTP(2); %%for bonn 0.16
        NC = FILTP(3); %2 for bonn 2
        STF1 = FILTP(4); %larger value smeers the time course!!! for bonn 0.57;
        STF2 = FILTP(4);
    end
end


filyert = 'fir';
%filyert = 'iir';

VF1ALL = [];
NVF1ALL = [];
xxf = 0; % for testing only!!!! should be zero!
ff = ff-xxf;
ff2 = ff2-xxf;


R = 'fir';

ssff = STF1;
AddFreq = 0;
%Filtype = 0; % for 1 'bandpass' for 0 use the new way 'filtfilt' changed on Aug 4 2021

Voriginal = V;
if removeHighLow % 2022 aug

    %remove artifact by filtering above ff2 and below ff1
    %than take the mean and subtract it from V
    AA = 40;
    ssfff = 0.98;
    DFF2 = 40;
    F = ff2+1*80;
    AddFreq = 0;
    global VF1UP VF1DOWN V_
    [VF1UP,DFilter] = bandpass(V,[AddFreq+F-DFF2 AddFreq+F+DFF2],1/dt,'ImpulseResponse',filyert, 'Steepness',ssfff,'StopbandAttenuation',AA);

    ssfff = STF1;
    AddFreq = 0;

    DFilterL = DFilter;
    DFF2 = 20;
    F = ff-2.0*DFF2;

    [VF1DOWN,DFilter] = bandpass(V,[AddFreq+F-DFF2 AddFreq+F+DFF2],1/dt,'ImpulseResponse',filyert, 'Steepness',ssfff,'StopbandAttenuation',AA);
    %factLow = %calculate the expected artifact at ff-2.5*DFF to that of ff
    %assuming a time constant of 10 ms without an electrode:
    Rati = abs((1+1i*ff*0.01)./(1+1i*(ff-2.5*DFF2)*0.01));
    %
    V_ = V-(1*VF1UP+1*VF1DOWN)/1;
    V = V_;
    % here we use VC from a previous run to subtract the V and see how good
    % it works.
    global VC;
    %V = V-VC/2;
end

A = 80;

if Filtype ==1;
    %% test for z outside the frequencies by 25H above. only for V.

    [VF1,DFilter] = bandpass(V,[AddFreq+ff-DFF AddFreq+ff+DFF],1/dt,'ImpulseResponse',filyert, 'Steepness',ssff,'StopbandAttenuation',A);

    IF1 = bandpass(I,[AddFreq+ff-DFF AddFreq+ff+DFF],1/dt,'ImpulseResponse',filyert, 'Steepness',ssff,'StopbandAttenuation',A);

    VF2 = bandpass(V,[AddFreq+ff2-DFF AddFreq+ff2+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF1,'StopbandAttenuation',A);

    IF2 = bandpass(I,[AddFreq+ff2-DFF AddFreq+ff2+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF2,'StopbandAttenuation',A);
    Filtmethod = 'bandpass...'
else
    firr = 1;
    if firr;

        DFF= 55;%44.25; %  2Hz %DFF 22.25 and 35 when filter z or g
        Windtype = 'ls';
        FilterOrderI = 5355;% 1500; %900 is good!!
        bpFilt_FF = designfilt('bandpassfir', 'FilterOrder', FilterOrderI,...
            'CutoffFrequency1', AddFreq+ff-DFF, 'CutoffFrequency2', AddFreq+ff+DFF,...
            'SampleRate', 1/dt);
        bpFilt_FF2 = designfilt('bandpassfir', 'FilterOrder', FilterOrderI,...
            'CutoffFrequency1', AddFreq+ff2-DFF, 'CutoffFrequency2', AddFreq+ff2+DFF,...
            'SampleRate', 1/dt);
        'hey...designfilt filtering ...'
    else %iir
        FilterOrderI = 100; %1000 is good!!
        bpFilt_FF = designfilt('bandpassiir', 'FilterOrder', FilterOrderI, ...
            'CutoffFrequency1', AddFreq+ff-DFF, 'CutoffFrequency2', AddFreq+ff+DFF,...
            'SampleRate', 1/dt);
        bpFilt_FF2 = designfilt('bandpassiir', 'FilterOrder', FilterOrderI, ...
            'CutoffFrequency1', AddFreq+ff2-DFF, 'CutoffFrequency2', AddFreq+ff2+DFF,...
            'SampleRate', 1/dt);
    end




    filtfiltt= 1;

    if filtfiltt ==1
        VF1 = filtfilt(bpFilt_FF,V);
        IF1 = filtfilt(bpFilt_FF,I);

        VF2 = filtfilt(bpFilt_FF2,V);
        IF2 = filtfilt(bpFilt_FF2,I);
        Filtmethod = 'filtfilt...'
    else
        VF1 = filter(bpFilt_FF,V);
        IF1 = filter(bpFilt_FF,I);

        VF2 = filter(bpFilt_FF2,V);
        IF2 = filter(bpFilt_FF2,I);

        Filtmethod = 'filter...'
    end

end


if Filtype == 3;
    fs = 1/dt;
    [wtt,fww] = cwt(V, fs);
    wdf = 8;
    VF1 = icwt(wtt,fww,[ff-wdf ff+wdf],'SignalMean',mean(V),'amor'); % bandpass between ff-2 to ff+2 Hz
    VF2 = icwt(wtt,fww,[ff2-wdf ff2+wdf],'SignalMean',mean(V),'amor'); % bandpass between ff-2 to ff+2 Hz

    VF1 = VF1-mean(VF1(round(searchtime(1)/dt):round(searchtime(2)/dt)));
    VF2 = VF2-mean(VF2(round(searchtime(1)/dt):round(searchtime(2)/dt)));
end

if VCfirst % processing VC
    V = VA; % 2022 Sep to remove VC from V before.
end

evf1 = envelope(VF1);
eif1 = envelope(IF1);

evf2 = envelope(VF2);
eif2 = envelope(IF2);

global DavidData;
if isempty(DavidData); %was not yet determined.
    DavidData = 0;
end

if DavidData
    '...david data'
    global VF1David VF2David comps_hatV comps_hatI V_or_I comps_hatV_IL Stable_demix_op
    fffknown = 0;% this is to use if we want to specify the frequencies as used (1) or used those measured above (0)
    if fffknown;
        ff = 145; ff2 =  260;
    end
    V_or_I = 1;
    %stable_demix
    
    demixmode = 'linear';
    demixmode = 'wt';
    %demixmode = 'wt2';
   Stable_demix_op  = 1;
    if Stable_demix_op;

        [comps_hatV,err] = stable_demix(V, 1/dt, 2*pi*ff, 2*pi*ff2,mode = demixmode);
    else
        [comps_hatV,err] = demix(V, 1/dt, 2*pi*ff, 2*pi*ff2);
    end

    VCdemix = comps_hatV(1,:);
    VF1 = comps_hatV(2,:);
    VF2 = comps_hatV(3,:);
    V_or_I = 0;
    if Stable_demix_op;
        [comps_hatI,err] = stable_demix(I, 1/dt, 2*pi*ff, 2*pi*ff2,mode = demixmode);
    else
        [comps_hatI,err] = demix(I, 1/dt, 2*pi*ff, 2*pi*ff2);
    end

    IC = comps_hatI(1,:);
    IF1 = comps_hatI(2,:);
    IF2 = comps_hatI(3,:);

end




hilbVf1 = hilbert(VF1);
hilbIf1 = hilbert(IF1);


hilbVf2 = hilbert(VF2);
hilbIf2 = hilbert(IF2);

meif2 = mean(eif2(2000:end-2000));
meif1 = mean(eif1(2000:end-2000));

if size(eif2)~= size(evf2);
    eif2 = eif2';
    hilbIf2 = hilbIf2.';
end
im2 = (evf2./(eif2)); % this is the impedance of the rs (at high freq), electrode resistance
%im2 = (evf2./(meif2));

global Zt;
if size(eif1)~= size(evf1);

    eif1 = eif1';
    hilbIf1 = hilbIf1.';
end
Zt = evf1./eif1; % total impedance including electrode for the low freq
%Zt = evf1./meif1; % total impedance including electrode for the low freq


%Rs_meas = im2; % vector of the Re based on high freq (aug 2020)
w = 2*pi*ff;

% if cleaning was done use the V from there, otherwise use use the original
V = Voriginal; % aug 7 2022 for the artifact removal



if modec ==1;  %calculate the c from phase analysis, change to take it from the higher freq
    % the assumption here that the for the higher frequency, at resting
    % condition g of the cell is much smaller than w*c and thus we solve the
    % problem as if there is only the resistor of the electtode and the
    % capcitor of the cell.
    % we assume: Z=~ Rs-j*(1/W*c) when approximating.
    % hence we will find the phase relationship between V and I and inverse
    % to find c.
    fcf = ff2;
    STF2 = 0.93;%0.6
    fi_2 = IF2;% bandpass(I,[fcf-DFF fcf+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF2);
    fv2 = VF2;% bandpass(V,[fcf-DFF fcf+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF2);
    fi_2 =  bandpass(I,[fcf-DFF fcf+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF2);
    fv2 =  bandpass(V,[fcf-DFF fcf+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF2);



    %fi = fi';

    s_fv = size(fv2);
    s_fi = size(fi_2);
    if s_fv ~= s_fi;
        fi_2 = fi_2';
    end

    fii2 = fi_2(round(searchtime(1)/dt):round(searchtime(2)/dt));
    fvv2 = fv2(round(searchtime(1)/dt):round(searchtime(2)/dt));

    global cmm2 ang11
    ang11 = median((angle(hilbert(fvv2))-angle(hilbert(fii2))));


    RRR = max(abs(fft(fvv2-mean(fvv2))))/max(abs(fft(fii2-mean(fii2)))); % this gives the best Jan 2022
    % RRR = max(real(fft(fvv2-mean(fvv2))))/max(real(fft(fii2-mean(fii2)))); % this should be correct (20/10/2021)
    %RRR = max(real(fft(V-mean(V))))/max(real(fft(I)-mean(I)));
    %RRR = max(real(fft(V-mean(V))))/max(real(fft(I-mean(I))));
    %2022   RRR = max(real(fft(fvv2-mean(fvv2))./fft(fii2-mean(fii2)))); % this should be correct (20/1/2022) % checked well!!!
    %comment although I thought that this (last line)is better it is not!!! I checked in simulations. Feb 24 2022
    cmm = c;
    factorC = 1;
    RRRR = RRR/1e6;
    %  cmm =  abs(1/(atan(ang11)*RRR*2*pi*ff2))*factorC; %with a mistake here
    cmm =  abs(1/(tan(ang11)*RRR*2*pi*ff2))*factorC; % this should be correct based on what Mor Ovadia found 20/10/2021)

    c = cmm*1.0;
    theta = ang11

end

%Gtotal = (w^2*c^2*(Zt.^2-Rs_meas.^2-1./(w^2*c^2)))./(2*Rs_meas); % transformation from Z to G

%exact


k = w*c;
Rsm = Rs_meas*1.0;
Zt = Zt*1;

NewExact3 = 1; %Dec 14 2020
if NewExact3% solving based on absolute value:
    %%Solve[Abs (r + 1/(g + I*w1*c)) == Abs (z1) && Abs (r + 1/(g + I*w2*c)) == Abs (z2), {r, g}]
    g1 = [];
    g2 = [];
    %test for imag part to be minus for current only.

    z1 = 1*hilbVf1./hilbIf1;
    z2 = 1*hilbVf2./hilbIf2;
    
    
    testiFFT =0;
    if testiFFT
        DFF = (ff2-ff)*0.11;
        atten = 1/30;
        DFFforVC = 70;

        [z1 z2 VC] = getZinverseFFT2022(V,I,dt,ff,ff2,4,DFF,atten,DFFforVC);
        4;
    end

    lowpassZ = 0;
    if lowpassZ;
        Filtype = 1;
        if Filtype == 0
            bpFilt_FFZ=  designfilt('lowpassiir', 'FilterOrder', 14, ...
                'PassbandFrequency', ff-40, 'PassbandRipple', 0.2,...
                'SampleRate', 1/dt);
            z1 = filtfilt(bpFilt_FFZ,z1);
            z2 = filtfilt(bpFilt_FFZ,z2);
        else
            z1 = lowpass(z1,ff,1/dt,'ImpulseResponse','iir','StopbandAttenuation',80,'Steepness', 0.999999);
            z2 = lowpass(z2,ff2,1/dt,'ImpulseResponse','iir','StopbandAttenuation',80,'Steepness', 0.999999);
        end

    end


    cleannoise = 0; % to smooth the
    if cleannoise;
        wws = 751;
        z11r = smooth(real(z1),wws);
        z11i = smooth(imag(z1),wws);
        z1o = z1;
        z1 = z11r+1i*z11i;
        z1 = z1';
        z22r = smooth(real(z2),wws);
        z22i = smooth(imag(z2),wws);
        z2 = z22r+1i*z22i;
        z2 = z2';


        %figure
        %plot(imag(z1)); hold on; plot(imag(z1o));
    end
    'done clean'

    Z1_orig = z1;
    Z2_orig = z2;
    if ~exist('BoostCe');
        boostcehere = 0;
    else;
        boostcehere = BoostCe;
    end

    st1 = round(searchtime(1)/dt);
    st2 = round(searchtime(2)/dt);
    zz1 = z1(st1:st2);
    zz2 = z2(st1:st2);
    global estce ce estce2
    estce2  = mean(abs((1./abs(zz1)-1./abs(zz2))./(2*pi*(ff-ff2))));
    estce  = mean(abs((1./(zz1)-1./(zz2))./(2*pi*(ff-ff2)))); %first the diff and then the abs


    if or(boostcehere,hyBridCe ==1); %correcting for electrode capac.
        gt1 = (1./z1);
        gt2 = (1./z2);
        % estce = ce;
        Kce = 1.0; %2021 should be 0.95
        estce = estce*Kce;
        gc11 = gt1-1i*estce*2*pi*ff;
        gc21 = gt2-1i*estce*2*pi*ff2;
        z1 = 1./gc11; % a new value for z1 and z2
        z2 = 1./gc21;
        %use this to find the correct c (cell cpacitance): Jan 11 2022,
        %checked and good!
        ang11 = mean(angle(z2(round(searchtime(1)/dt):round(searchtime(2)/dt))));
        RRR = mean(real(z2(round(searchtime(1)/dt):round(searchtime(2)/dt))));
        cmm =  abs(1/(tan(ang11)*RRR*2*pi*ff2))*factorC;
        c = cmm;
    end

    w1 = 2*pi*ff;
    w2 = 2*pi*ff2;
    z1a = (z1);
    z2a = (z2);
    if hyBridCe; % here we use the value of z before correction with pipette capac, but use later c as was found with pippete correction
        z1a =  Z1_orig;
        z2a =  Z2_orig;
    end

    %%%
    %% this comes from solving:
    %

    if compensLowRin % try sep 2022;
        z1aaa = z1a;
        z2aaa = z2a;

        if 1%resale only when z2>z1;
            crossingindices = find(real(z2a)>real(z1a));
            bymuch = real(z2a(crossingindices))-real(z1a(crossingindices));
            z2areal = real(z2a);
            z2aimag = imag(z2a);
            z2areal(crossingindices)= z2areal(crossingindices)-bymuch*1.1;
            z2a =z2areal+1i*z2aimag;
            z1 = z1a;
            z2 = z2a;
        end
        if 0
            z1aa = z1a(5000:end-5000);
            z2aa = z2a(5000:end-5000);
            moveby = max(real(z2aa)-real(z1aa));
            if moveby>0;
                z2areal = real(z2a);

                z2areal = z2areal-moveby*1.5;;
                z2aimag = imag(z2a);
                z2a =z2areal+1i*z2aimag;
                z1 = z1a;
                z2 = z2a;
            end
        end
    end
    rcc = 0; %try to correct jan 2023
    if rcc
    z1a = 1*z1a+30*MC;
    z2a = 1*z2a-30*MC;

    end
    Rs_meas = (1/(2*(1i*c*w1 - 1i*c*w2)))*(1i*c*w1*z1a - 1i*c*w2*z1a + 1i*c*w1*z2a -...
        1i*c*w2*z2a + ((-1i*c*w1*z1a + 1i*c*w2*z1a - 1i*c*w1*z2a + 1i*c*w2*z2a).^2 -...
        4*(1i*c*w1 - 1i*c*w2)*(z1a - z2a + 1i*c*w1*z1a.*z2a - 1i*c*w2*z1a.*z2a)).^0.5);
    if rcc  %jan 2023
        
    end

    Gtotal = 1*((1i*(1i + c*Rs_meas*w1 - c*w1*z1a))./(Rs_meas - z1a));
    
    %setting the z1 and z2 as the last ones. 2023
    z1 = z1a;
    z2 = z2a;

    newsol2022 = 1;
    if newsol2022;% found today (7/4/2022) that mathematica (application) gives another way:

    end

    global GTTT;
    GTTT = Gtotal;


    lowpassG = 0; % gal check this part

    if lowpassG;
        if Filtype == 0 %when using 'filt' to filter V
            bpFilt_FFG=  designfilt('lowpassiir', 'FilterOrder', 6, ...
                'PassbandFrequency', ff-20, 'PassbandRipple', 0.8,...
                'SampleRate', 1/dt);
            Gtotal = filtfilt(bpFilt_FFG,Gtotal);
            Rs_meas = filtfilt(bpFilt_FFG,Rs_meas);
        else
            Gtotal = lowpass(Gtotal,ff-20,1/dt,'ImpulseResponse','iir','StopbandAttenuation',80,'Steepness', 0.999999);
            Rs_meas = lowpass(Rs_meas,ff-20,1/dt,'ImpulseResponse','iir','StopbandAttenuation',80,'Steepness', 0.999999);
        end

    end


end
re = Rs_meas;

%%% Jan16_2021 boost for Ce electrode:
% start by trying to add to  changes in re due to ce:
% then compute again z1 and z2


%reconstruct Ge and Gi
DFF = 5;
%cleaning the voltage and the total conductance.
global VC VCC1
VCC1 = VC;
if 1 % processing VC
    filyert = 'fir';
    DFF = 60;

    VC = V;
    Iclean = I;
    LCC = lc;
    A = 25;
    if 1 % essential do not change

        for MF = 1:length(lc)/2;
            ff3toend = lc(MF)*df+99;
            STF3 = 0.54;
            VC = bandstop(VC,[ff3toend-DFF ff3toend+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',A);
            Iclean = bandstop(Iclean,[ff3toend-DFF ff3toend+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',A);
            global VC1
            VC1 = VC;
        end
    end

    %lowpass the VC and Iclean below the lower frequency of the two
    %frequencues Jan 2022;

    if 0
        if Filtype == 1 %when using 'filt' to filter V
            VC =  lowpass(VC,ff-10,1/dt,'ImpulseResponse',filyert,'Steepness',0.997,'StopbandAttenuation',A);
            Iclean =  lowpass(Iclean,ff-10,1/dt,'ImpulseResponse',filyert,'Steepness',0.9997,'StopbandAttenuation',A) ;
            VC =  lowpass(V,ff-30,1/dt,'ImpulseResponse','iir','Steepness',0.999,'StopbandAttenuation',120);
            Iclean =  lowpass(I,ff-30,1/dt,'ImpulseResponse','iir','Steepness',0.999,'StopbandAttenuation',120) ;
            global VC2
            VC2 = VC;

        else % using the 'filt' option
            bpFilt_FFVclean=  designfilt('lowpassiir', 'FilterOrder', 10, ...
                'PassbandFrequency', ff-10, 'PassbandRipple', 0.2,...
                'SampleRate', 1/dt);
            VC = filtfilt(bpFilt_FFVclean,VC);
            Iclean = filtfilt(bpFilt_FFVclean,Iclean);
            'lastone....'
        end

    end

    global ICLEAN;
    ICLEAN = Iclean;

    if 1 % changed to zero 14/4/2021 from 1
        Gtotal = real(Gtotal);
        B= 100;
        %Gtotal = bandstop(Gtotal,[ff-DFF ff+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',B);
        %Gtotal = bandstop(Gtotal,[ff2-DFF ff2+DFF],1/dt,'ImpulseResponse',filyert,  'Steepness',STF3,'StopbandAttenuation',B);
        %        Gtotal = lowpass(Gtotal,ff-60,1/dt,'ImpulseResponse',filyert,'Steepness',0.993,'StopbandAttenuation',B);
    end
end
sVV = size(V);
sVF1 = size(VF1);
size(VF2);

if sVV ~= sVF1;
    V = V';
end
sizVC = size(VC);
%VC = V-(VF1+VF2); %2021 may check!

if DavidData %using the method of David Ungarish
    VC = VCdemix;
end
% dvdt and c
if 1 % 2022 to xxxx checkit !!!!
    svv = size(VC);
    sii = size(Iclean);
    if svv~=sii;
        Iclean = Iclean';
    end

    VC = VC-Iclean.*abs(re); % to bridge balance the voltage
    VC = VC-Iclean.*real(re); % to bridge balance the voltage %changed on 4/11/2021
end
sssvc = size(VC)
sdv = size(VC);
if sdv(1) >1;
    VC = VC';
end
svc = size(VC);

dv = diff(VC);
global DVV;

sdv = size(dv);
if sdv(1) >1;
    dv = dv';

end
DVV = dv;
cdvdt = c*[dv dv(end)]./dt;
%% test for problems with VC using the simulations only
testSIM = 0; % we use the VC and the real GT to see how much the filtering of the voltage to create VC affect the measurement
if testSIM
    global GTT;
    Gtotal = GTT; %GE+GI from the main simulator
    % running this show that the filters that generated GT are not doing much
    % artifact.
end

%finding the resting potential and gleak:
t1points = round(searchtime(1)/dt);
t2points =round(searchtime(2)/dt);


Vwind = VC(t1points:t2points);
%Vwind = Gtotal(t1points:t2points); %aug 2022 with Sagi!!!
% because of the oscillations in G, this last change can be a problem...
%so I took it back.
l5p = prctile(Vwind,15); % 5%
l5p = prctile(Vwind,15); % 5%%aug 2022 changed with Sagi!!!


viniceslow5 = find(Vwind<=l5p)+t1points;
global VRR;

vr = mean(VC(viniceslow5));
g0 = mean(Gtotal(viniceslow5));

%caluclating the ge and gi from clean v and from total g as in lampl
%2019 first figure and equations.
global REVS
REVS = reves;
gsyn = Gtotal-g0;
global GSYN
GSYN = gsyn;
GT = Gtotal;
ve = reves(1);
vi = reves(3); %confusing that reves should have 3 entries (~DANIEL)
vl = vr; %which is also vleak.
%vll= VC-vl;
vll= VC-(vr);
vee = VC-ve;
vii = VC-vi;
TT2 = vii;


dvdts = size(cdvdt);
if size(gsyn) ~= dvdts;
    gsyn = gsyn';
end

svee  =size(vee);
svii = size(vii);
%g0 = g0.*(gsyn/(g0));
%g0 = gsyn.*g0./(g0+gsyn);
if ~exist('cableBoost');
    boostcable  = 0;
else
    boostcable = cableBoost;
end

if boostcable % this is boosting for cable case.
    global Gnew
    g0old = g0;
    %kg = mean(gsyn/g0);
    kg = abs(mean(gsyn(15000:end-15000)/g0));
    g0 = g0old*(1-exp(-1*(gsyn./g0).^2)); %checked for change from -4 to -2 aug 4 2021
    g0(g0<0) = 0;
    'cableboosting....'

    Gnew = g0;
end

if size(Iclean) ~= size(VC);
    Iclean = Iclean';
end

%gsyn = Gtotal'-g0;%% test 2022 2_13_2022 to see if we can corrrect the synap conductance after correcting gleak
gi = (1.0*cdvdt+1*g0.*vll+gsyn.*vee-1*Iclean)./(vee-vii); %inhibition
gl = g0;
if boostcable;
    gl = g0old;
end
ge = gsyn-gi;
ge = real(ge); % 2022 orig real
gi = real(gi);


Iin = gi.*vii;
Iex = ge.*vee;

lowpassing = 0;
if lowpassing;
    ge = lowpass(ge,12,1/dt,'Steepness',0.9999999,'StopbandAttenuation',120);
    gi= lowpass(gi,12,1/dt,'Steepness',0.9999999,'StopbandAttenuation',120);
    GT = ge+gi+g0;
end


if plotit
    figure
    pointst = 1;
    sp1 = subplot(4,1,1);
    Vs = V(pointst:end-pointst);
    tt = pointst*dt;
    plot([1:1:length(Vs)]*dt,Vs); hold on; % plot([1-tt 2-tt],[-0.04 -0.04],'b')
    ylabel('recorded Vm (V)')
    xlabel('Time (S)')
    hold on;
    sp2 = subplot(4,1,2);
    gee = sgolayfilt(ge(pointst:end-pointst),1,15);
    gii = sgolayfilt(gi(pointst:end-pointst),1,15);
    VCC = sgolayfilt(VC(pointst:end-pointst),1,15);
    env = sgolayfilt(evf1(pointst:end-pointst),1,15);
    env2 = sgolayfilt(evf2(pointst:end-pointst),1,3);

    env = (evf1(pointst:end-pointst));
    env2 = (evf2(pointst:end-pointst));

    env21 = env-1*env2;
    %env21 = lowpass(env21,3,1/dt,'Steepness',0.999);
    if 1
        % gee = bandstop(gee,[ff-DFF ff+DFF],1/dt);
        % gee = bandstop(gee,[ff2-DFF ff2+DFF],1/dt);

        %gii = bandstop(gii,[ff-DFF ff+DFF],1/dt);
        %gii = bandstop(gii,[ff2-DFF ff2+DFF],1/dt);
        timestartend = 0.15;
        gee(1:round(timestartend/dt)) = 0;
        gee(end-round(timestartend/dt):end) = 0;

        gii(1:round(timestartend/dt)) = 0;
        gii(end-round(timestartend/dt):end) = 0;
        plot([1:1:length(gee)]*dt,gee,'g'); hold on;  plot([1:1:length(gee)]*dt,gii,'r');
        hold on;
        plot([1:1:length(gee)]*dt,gii+gee,'k');
        ylabel('Cond (S)')
        xlabel('Time (S)')

    end
    hold on;

    sp3 = subplot(4,1,3);
    plot([1:1:length(VCC)]*dt,VCC,'k');

    ylabel('Vm (V)')
    title('Filt Vm');
    xlabel('Time (S)')
    % plot(GE(1000:end-1000),'g'); hold on; plot(GI(1000:end-1000),'r')


    sp4 = subplot(4,1,4);
    env = envelope(abs(z1));
    env2 = envelope(abs(z2));
    plot([1:1:length(env)]*dt,env,'k');
    hold on;
    plot([1:1:length(env2)]*dt,env2,'k');
    xlabel('Time (S)')

    linkaxes([sp1 sp2 sp3 sp4],'x');

    figure
    hold on;
    MC = 1000000;
    MC = 1;
    %plot([1:1:length(evf1)]*dt,real(z1)./1e6,'k');

    plot([1:1:length(re)]*dt,abs(re)./1e6,'b');
    plot([1:1:length(GT)]*dt,real(1./GT)./1e6,'m');
    plot([1:1:length(GT)]*dt,abs(1./GT)./1e6,'k');
    legend({'Re','realRin (cell)','absRin (cell)'})
    ylabel('Impedance inMohm')
    xlabel('Time (S)')
    set(gca,'ylim',[0 1000]);
    title('New')
    % plot([1:1:length(env)]*dt,MC*env21,'b');
end
