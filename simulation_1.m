function ret = simulation_1(stm_E)

arguments
    stm_E = StmInfo();
end

SPT = [.2:0.25:2, 2.8, 3.2];

stm_I = stm_E;
stm_I.trec = 4;


Kg = 1500;

sfdt = 30000; %40000 Hz %main dt setting
dt = 1/sfdt; %simulation time constant
VVV = []; %check noise 2023
ratiomeasured = [];
Freqtests = [];
Retests = [];
%create a slow varing changes in electrode resistance:
NoiseR = randn(1,500000);
NoiseR =  smooth(NoiseR,3311);
NoiseR = NoiseR(1000:end-1000);
NoiseR = sgolayfilt(NoiseR,1,3311);
NoiseR= NoiseR-mean(NoiseR);
NoiseR = NoiseR./std(NoiseR);
NoiseV = randn(1,500000);
NoiseV = sgolayfilt(NoiseV,1,15);
RRR = [];

MC = 1000000; % this is a scaling factor for moving from Mohm to Ohms.
t = [0:0.0001:2];
c = 1*1.5*1e-10;%  1*1.5*1e-10; ?????
tau = 0.0157;%0.01 %this is the SYNAPSE tau!
gt = t.*exp(-t/tau);
%taum  = 0.009; %tau of membrane 10 mstaum
Rcell = 30; %150 %2022 when Rin is too small ??????? ????? ?? ??? z1 and z2 intersect and this leads to problems.
gl = 1/(Rcell*MC); %gleak for 100MOhm input resistance
GLL = gl;
% Rinn = (1/gl)/MC;
GLL = gl;
taum = Rcell*MC*c;
Ge = 250*gl;
Gi = 460*gl;
vl = -0.07;
ve = -0.00;
vi = -0.08;
global GL
GL = gl;

for mm = 1:3

    V = [];
    VM = [];
    GE = [];
    GI = [];
    GHHN = [];
    GHHK = [];
    SII1 = [];
    SII2 = [];
    IClamp = [];
    ii = 0;

    Vp = vl;
    v = vl;

    ti = 0;
    te = 0;
    T = [];
    Ig = 0;


    taue = 0.001;
    gelec = (1/3)*gl;
    re = 0.05*1/gelec;
    ce = taue*(1/re);




    %     u = 0.7;
    %     trec = 0.6;
    %     taui = 0.1;
    %     ase = 1;
    %     freq = 50;
    %     p = 1;
    %     U = p*u;

    %A = [SPT 2.8 3.2];
    %SPT = A;


    epscampsE = ShortTerm_MTmodel(SPT,stm_E.u,stm_E.trec,stm_E.taui,stm_E.ase);
    epscampsI = ShortTerm_MTmodel(SPT,stm_I.u,stm_I.trec,stm_I.taui,stm_I.ase);

    if length(epscampsI) >= 11
        % LEGACY. not sure why we're doing this..
        epscampsI(9) = 0;
    end

    IIm = [];
    TTT = [];
    stimnumP= -1;
    %%%% some of the variable below are not used anymore as they were used
    %%%% for VC.
    k = 0;
    ii = 0;
    I = 0;
    Ic = 0;
    Im = 0;
    vprev = vl;
    vprevS = vl;
    IClampe = [];
    IClampi = [];

    IClampP = [];
    IClampL = [];
    perrorP = 0;
    perrorL = 0;
    perror = 0;
    Icp  = 0;
    Icl = 0;
    vm = vl;
    dv = 0;
    dvm = 0;
    Ice = 0;
    Ici = 0;
    Ig = 0;
    Ige = 0;
    Igi = 0;
    kkk = 0;
    Igl = 0;
    Igp = 0;
    l = 0;
    ALLC= [];
    III = [];
    vP = 0;
    dv = 0;
    Zexp = [];
    Iinj = [];

    VCevery = 0.0001; %sample and clamp every 1 ms..0.0004.
    counttimeVCevery = 0;
    %tt = 0:dt:4;
    giii = 0*gl/2;
    geee = 0*gl/8;
    TTT = 0:dt:1.5;
    % prenoise = randn(1,length(TTT))*0.0013;
    for tt = 0:dt:1.5;% 0:dt:4; ??? ???? this is the main time
        gi = giii;
        ge = geee;
        ii = ii+1;
        if mod(tt,0.01)== 0;
            ttt= tt;
        end

        stimnum = find(SPT(1:end-1)<=tt & SPT(2:end)>tt);
        if ~isempty(stimnum);
            ampe = epscampsE(stimnum);
            ampi = epscampsI(stimnum);

            if stimnum~=stimnumP
                te = 0;
                ti = 0;
                k = 0;
            end
            k = k+1;
            te = te+dt;
            if k <= round(0.003/dt); % 4 ms delay of I compared to E
                ti = 0;
                ampi = 0;
            end

            ti = ti+dt;
            if ti>0;

            end
            facttime = 1;
            factamp = 1;
            ge = factamp*1*8000*ampe*Ge*te^3*exp(-te/(tau*facttime))+geee;
            gep = ge;
            gi = factamp*0.4*8000*ampi*Gi*ti^3*exp(-ti/(tau*facttime))+giii;

            if t>2.1 %test Jan 2020
                gi = 0;
            end
            gip = gi;
            %trise = 0.009;

            %tdecay = 0.01;
            %ge = -ampe*Ge*1*((exp(-te/trise)-exp(-te/tdecay)));
            %gi = -ampe*Gi*1*((exp(-ti/trise)-exp(-ti/tdecay)));

            stimnumP = stimnum;
        else;
            ge = geee;
            gi = giii;
        end

        %-------------

        if ii> round(0.8/dt) & ii<= round(0.9/dt); %2.5 to 3.8 2022 and 3.6
            %if ii> round(3.0/dt) & ii<= round(3.4/dt);
            gi= 0*7*gl*0.5;%+1*1*giii+gip;
            ge = 0.75*gl*0.5;%+1*geee+gep;
            gi= ge*2; %0*7*gl*0.5+1*1*giii+gip;

            %ge = geee;
        end

        I = 0;
        SW = 1; %sine wave injection is on
        if SW == 1;
            Fss = [121.0 263 251 301 4411]; % 220 277
            Fss = [111 222*1 251 301 4411]*4; % 220 277 Feb2023
            %Fss = [103 196 251 301 4411]; % 220 277
            fsin1 = Fss(1);
            fsin2 = Fss(2);
            fsin3 = Fss(3);
            fsin4 = Fss(4);
            fsin5 = Fss(5);
            if mm == 1;
                ampsinfig = 0;
            else
                ampsinfig= 1;
            end

            ampSin1 = 1*ampsinfig*3*250*1*1e-12;
            SI1 = 1*ampSin1*sin(tt*2*pi*fsin1);
            SI2 = 1*ampSin1*sin(tt*2*pi*fsin2+0);
            SI3 = 0.0*ampSin1*sin(tt*2*pi*fsin3);
            SI4 = 0.0*ampSin1*sin(tt*2*pi*fsin4);
            SI5 = 0.0*sin(tt*2*pi*fsin5);
            SII1(ii) = SI1;
            SII2(ii) = SI2;
            % SI = SI1+SI2+SI3+SI4+SI5;
            %SI = SI1;
            if 1
                if ii>round(3.7/dt);
                    SI = SI2;

                end

            end
            SI = SI1+SI2;%+SI2+SI3;
            SI = SI1+SI2+SI3+SI4;
        else
            SI = 0;
        end

        vP = v;
        bw = 1;

        Ielec = 0*1*ce*dv/dt;
        %Im = I-Ielec+SI-0;
        Ipulse=0;
        Im = I+SI+Ipulse+Ielec;

        %Im = I+SI;
        %vm = v-Im*re;
        kkc = 1;
        scfacorG = 3;  %3 is in the figures. %feb2023

        kge = 2;%3 for paper
        kgi = 3; %3
        %vi = vl;

        %re = 1*(Rcell/4)*MC+rand*0*MC;

        %re = 0;


        %c = CCc*[1+gi*20];Re
        %c = 0.21*gi+CCc;
        CCC(ii) = c;
        CN = c;
        %Im = 0;
        % vm = vm+1*dvm;
        %% Sep 2020 add one more resistor Rp and Ce for the Rp capcaitance.
        %Vp is the voltage above the patch electrode, and v is above the
        %device resistor.
        EC = 0; %with (1) or without electrode capacitance (0).
        Rp = 1*EC*16*1.00000*MC; %this is the resistance directly connected to the cell. The capcitance goes above it.
        %Rp cannot be zero when EC = 1
        %tauelect = 0.0001;
        cenew = 0.6e-11; % for sylgard 3e-12 and regular 3e-11
        if mm == 2;
            ampre = 0;
        else
            ampre = 1;
        end

        re = ~EC*ampre*(30.00)*MC+rand*0*MC; % 30 for the paper
        RE = re;
        %% add noise to re:
        noiseRe = 0;
        if noiseRe == 1;
            re = re+ NoiseR(ii)*0.3*MC;
        end

        ce = cenew;
        %ce = tauelect/Rp;
        Taucell = c/gl;
        Taue = ce*re;
        Im = Im + 0*250*1*1e-12;
        if EC>0;
            if Rp>0
                Ic = Im-(Vp-vm)/Rp; % Im-(Vp-vm)/Rp;
            end
        else
            Ic = 0;
        end
        vm = vm+0*randn*0.00005+0*NoiseV(ii)*0.00006;
        VVV(ii) = vm;
        %  gi = 0;
        %adding HH conductances;
        HH  = 0;
        VIh = -0.045;
        GNAHH = 0;
        GKHH = 0;
        scaleGIh  = 1;

        GNAHH = 0;
        GKHH = 0;
        goutIh = 0;
        VNa = 0.05;
        VK = -0.09;



        % ve = vl; %% test only
        %vi = vl; %% test only
        %dvm = -dt*(1/(kkc*c))*(gl*(vm-vl)+1*kge*scfacorG*ge*(vm-ve)+1*kgi*scfacorG*gi*(vm-vi)-(Im-EC*Ic));
        dvm = -dt*(1/(kkc*c))*(gl*(vm-vl)+1*kge*scfacorG*ge*(vm-ve)+1*kgi*scfacorG*gi*(vm-vi)+0*GNAHH*(vm-VNa)+0*GKHH*(vm-VK)+0*scaleGIh*goutIh*(vm-VIh)-(Im-EC*Ic));
        dVp = EC*Ic*(dt/(ce));
        allc = gl*(vm-vl)+ge*(vm-ve)+gi*(vm-vi);

        %%

        vm = dvm+vm;
        v = vm;

        if EC == 0 %if EC = 0 no capa
            Vp = vm+Im*Rp;
        else
            Vp = Vp+dVp;
        end

        RR = re;
        RRR(ii) = re;
        v = 1*Vp+1*Im*RR;


        %%
        %%%
        % allc = ge*(vm-ve)+gi*(vm-vi);
        %compute theor:
        CC = c;
        gz = (gl+(kge*ge+kgi*gi)*scfacorG);
        bz = gz^2+(2*pi*Fss(1)*CC)^2;
        %Zexp(ii) = (re^2+(2*re*gz)/bz+1/bz)^0.5;
        wwc = 2*pi*Fss(1)*CC;
        Zexp(ii) = abs(re +gz/(gz^2+wwc^2)-i*wwc/((gz^2+wwc^2)))./MC;
        ZZ = Zexp;

        %dv = v-vP;
        IIm(ii)= Im;
        ALLC(ii) = allc;
        %%%%%%%%%%%%
        Iinj(ii) = Im; %changed from SI Dec 2020

        %voltage clamping:
        vprev = v;
        V(ii) = v;
        VM(ii) = vm;
        GE(ii)= kge*scfacorG*ge;
        GI(ii) = kgi*scfacorG*gi;

        T(ii)= tt;
        TTT(ii) = te;
        IClamp(ii) = Ic;
        IClampP(ii) = Icp;
        IClampL(ii) = Icl;

        IClampPS(ii) = Icp;
        IClampLS(ii) = Icl;


    end
    if abs(ampSin1) == 0
        VnoSin = V;
    end

    if mm == 2; %abs(re) == 0
        Vnoe = V;
    end


end



FILTP = [25 0.0001 3 0.79];

%FILTP = [15 0.0001 3 0.73]; %J?? 2022
Iinj = circshift(Iinj,0);

revs = [ve vl vi]; % [ve vl vi]
% to see the fields of the new function:
%function [ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2] = findGeGesimpleNew_cyclingFull4components_newBoostCe(V,I,sf,c,reves,searchtime,plotit,FILTP,cValue,BoostCe,cableBoost);
shifttt = 0*round(-0.00016/dt);
IinjTEMP = Iinj;
Iinj = circshift(Iinj,shifttt);

GE2 = [];
GI2 = [];
FrNum = 2;
combsF = FrNum*(FrNum-1)/2;
FK = 0;

%sooth V
Vor = V;
Iijnor = Iinj;
if 0
    if 0
        osg = 3;
        fsg = 217;
        Vss = sgolayfilt(V,osg,fsg);
        Iinjss = sgolayfilt(Iinj,osg,fsg);
    else
        nn = 11;
        deg = 5;
        Vss= smooth(V,nn,'sgolay',deg) ;
        Iinjss = smooth(Iinj,nn,'sgolay',deg) ;
    end
    V = Vss;
    Iinj = Iinjss;


else
    V = Vor;
    Iinj = Iijnor;

end

lowpasss = 0;
if lowpasss
    Vss= lowpass(V,Fss(2)+60,1/dt,'Steepness',0.97);
    V = Vss;

end

ret = struct();
ret.V = V;
ret.I = Iinj;
ret.t = T;
ret.gt = struct(GE=GE, GI=GI, GL=GL, GTOT=GE+GI+GL);

ret.cell_info = cellInfo(c=c,ve=revs(1),vi=revs(3));
ret.expr_info = exprInfo(Fs=1/dt,ws=2*pi*[fsin1,fsin2]);


% -------------
