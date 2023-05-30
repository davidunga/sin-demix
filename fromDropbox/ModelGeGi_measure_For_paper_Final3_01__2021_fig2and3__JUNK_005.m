%Simulating a cell with electrode in series. The cell recive depressing
%excitatory and inhibitory synaptic inputs and using AC analysis we reveal
%the changes in total condctance. The program is also for close loop (not
%active).
%units are in nA, ms, mV, Mohm
%close all;
global DavidData %Feb 2023
DavidData = 1; %feb2023
global Kg %scale used in fitting functios later
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
for HHH = 1:1;
    %close all
    
    figure
    for mm = 1:3;
        MC = 1000000; % this is a scaling factor for moving from Mohm to Ohms.
        t = [0:0.0001:2];%4
        %%%%%%%%%%%%
        c = 1*1.5*1e-10;%  1*1.5*1e-10; ?????
        %%%%%%%%%%%%
        tau = 0.0157;%0.01 %this is the SYNAPSE tau!
        gt = t.*exp(-t/tau);
        
        %taum  = 0.009; %tau of membrane 10 mstaum
        Rcell = 30; %150 %2022 when Rin is too small ??????? ????? ?? ??? z1 and z2 intersect and this leads to problems. 
        gl = 1/(Rcell*MC); %gleak for 100MOhm input resistance
        GLL = gl;
        % Rinn = (1/gl)/MC;
        GLL = gl;
        taum = Rcell*MC*c
        % c = taum/(Rcell*MC); %capac of the cell
        
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
        Ge = 250*gl;
        Gi = 460*gl;
        % ge = gl*3;
        % gi = gl*5;
        vl = -0.07%-0.07
        ve = -0.00; %0
        vi = -0.08; %-0.08
        Vp = vl;
        %ve = vl;
        %vi = vl;
        v = vl;
        global GL
        GL = gl;
        
        ti = 0;
        te = 0;
        T = [];
        Ig = 0;
        
        
        taue = 0.001;
        gelec = (1/3)*gl;
        re = 0.05*1/gelec;
        ce = taue*(1/re);
        %ce = 0;
        %  re = 0*(1/gl)/122; %will change later...
        %----------
        
        
        u = 0.7;
        trec = 0.6;
        taui = 0.1;
        ase = 1;
        freq = 50;
        p = 1;
        U = p*u;
        SPT = [ 1:0.125:2];
        A = [SPT 2.8 3.2];
        SPT = A;
        
        [epscampsE SPT] = ShortTerm_MTmodel(SPT,u,trec,taui,ase);
        u = 0.7;
        trec = 4;
        taui = 0.1;
        
        [epscampsI s] = ShortTerm_MTmodel(SPT,u,trec,taui,ase);
        
        epscampsI(9) = 0;
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
        
        IClampP = []
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
            if 0
                if ii> round(3.4/dt) & ii<= round(3.8/dt);
                    gi= 0*gl*0.5;
                    ge =0*gl*0.5;
                end
            end
            
            if ii> round(0.8/dt) & ii<= round(0.9/dt); %2.5 to 3.8 2022 and 3.6
                %if ii> round(3.0/dt) & ii<= round(3.4/dt);
                gi= 0*7*gl*0.5;%+1*1*giii+gip;
                ge = 0.75*gl*0.5;%+1*geee+gep;
                gi= ge*2; %0*7*gl*0.5+1*1*giii+gip;

                %ge = geee;
            end
            
            Ipulse =0;
            if Ipulse
                if (ii >round(0.7/dt)) & (ii <= round(3.6/dt));
                    Ipulse = -3*50*1e-12;
                else
                    Ipulse = 0;
                end
            end
            
            
            
            
            I = 0*Ipulse;
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
            if 0
                'ffffffffffffffffffff'
                v = re*Im + vm +0.000*rand ;
                vm = vm+h*dvm;
            end
            
            Ielec = 0*1*ce*dv/dt;
            %Im = I-Ielec+SI-0;
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
            if HH
                
                if ii == 1;
                    nhh = 0;
                    mhh = 0;
                    hhh = 0;
                    mhii = 0;
                end
                %[goutNa goutK nhh mhh hhh] = getNa_K_HH(ii,vm*1000,dt,nhh,mhh,hhh);
                [goutIh  mhii] = get_Ih_Current(ii,vm,dt,mhii);
                scaleGIh = 1/550000;
                %GNAHH = goutNa*gl*40;
                %GKHH = goutK*gl*5;
                %VNa = 0.05;
                
                %VK = -0.09;
                %GHHN(ii) = GNAHH;
                %GHHK(ii) = GKHH;
                GHHK(ii) = goutIh;
                %             if ii ==600;
                %                 return
                %             end
                
            else
                GNAHH = 0;
                GKHH = 0;
                goutIh = 0;
                VNa = 0.05;
                VK = -0.09;
            end
            
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
            
            if 0 %test for the effect of stimulus artifact
                if tt > 0.98 & tt<=0.981
                    if tt == 0.98
                        vvv = v;
                    end
                    v = v+0.001;
                    if tt == 0.981;
                        v = vvv;
                    end
                end
            end
            
            
            
            
            
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
end

if 0;%add somenoise to V %David Feb2023
V = V+prenoise;
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



if FrNum>2;
    for Fii = 1:combsF;
        for Fjj = Fii+1:FrNum;
            FreqArray= [Fii Fjj];
            FK = FK+1;
            FiltType = 0; %1 for bandpass, 0 for filtfilt
            [ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGi_MultiFreq_v003(V,Iinj,1/dt,c,revs,[0.3 0.5],0,FILTP,-1,0,0,FreqArray,FiltType);
            GE2(FK,:)= ge;
            GI2(FK,:)= gi;
        end
    end
    ge = mean(GE2);
    gi = mean(GI2);
else
    FreqArray= [1 2];
    FiltType = 0; %1 for bandpass, 0 for filtfilt 3 for wavelet
   % [ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGi_MultiFreq_v002(V,Iinj,1/dt,c,revs,[0.3 0.47],0,FILTP,-1,0,0,FreqArray);
    hybridCe = 0; 
    BoostCe = 0;
    %V = smooth(V,'sgolay');
    removeHighLow = 0; %aug 2022 new way to remove artifact
    FILTP = [30 0.0001 3 0.83];
    compensLowRin = 0; %to prevent z2>z1, lowering the real z2 
    %try to resample to higher SF and then resample back after. 
    if 0
    RN = 20;
    METHOD= 'pchip';
    Tx = [1:1:length(T)*RN]*(dt/RN);
    VR = resample(V,T,RN/dt,1,RN,METHOD);
    IinjR= resample(Iinj,T,RN/dt,1,RN,METHOD);
    GER = resample(GE,T,RN/dt,1,RN,METHOD);
    GIR = resample(GI,T,RN/dt,1,RN,METHOD);
    dtR = dt/RN;

    [geR,giR,glR,reR,VCR,GTR,ZtR,cmmR,XfoundR,ff,ff2,g1R,g2R,z1R,z2R,vl,Iin,Iex]=  findGeGi_MultiFreq_v005_temp(VR,IinjR,1/dtR,c,revs,[0.4 0.7],1,FILTP,-1,BoostCe,0,FreqArray,FiltType,hybridCe,removeHighLow,compensLowRin);% dsds ?????%     
    %
    
    else

    cell_info = cellInfo(c=c,ve=revs(1),vi=revs(3));
    expr_info = exprInfo(Fs=1/dt,ws=2*pi*[fsin1,fsin2]);
    res = estimate_conductances(Iinj,V,expr_info,cell_info,stable_demix=false);
    return

    [ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGi_MultiFreq_v005_temp(V,Iinj,1/dt,c,revs,[0.4 0.7],1,FILTP,-1,BoostCe,0,FreqArray,FiltType,hybridCe,removeHighLow,compensLowRin);% dsds ?????%     
  % [ge3,gi3,gl3,re3,VC3,GT3,Zt3,cmm3,Xfound3,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGi_MultiFreq_v003(V,Iinj,1/dt,c,revs,[0.4 0.7],1,FILTP,-1,BoostCe,0,FreqArray,FiltType);
    %
    end

    
end
%[ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGesimpleNew_cyclingFull4components_newBoostCe(V,Iinj,1/dt,c,revs,[0.3 0.47],0,FILTP,-1,0,0);

return
GTOT = GT;
%return
%findGeGesimpleNew_cyclingFull4components_newBoostCe
global evf1 evf2 VF1 VF2 IF1 IF2 eif1 eif2
%ploting panels of figures:


coV = 1000;
coI = 1e9;
coG = 1e9;
coZ = 1/MC;
%%% voltage real
FS = 7;
fontA = 10;
p1a = 0.19;% 0.17
p2a = 0.19;% 0.12

colorz1 = [0 0.4470 0.7410];
colorz2 = [0.8500 0.3250 0.0980];
if 0
    % two stages: take the coordinates of the destination, copy into f2 the
    %original figure and then set to the correct position.
    %ax1pos = get(sp1,'position');% from the multipanel figure
    %newAxcopied = copyobj(ax1,f2)
    %set(newAxcopied,'position',ax1pos) %to the
    
end



fbig = open('skel_figure2_method2_nonfilled_f2part1.fig');


figure(110)
ax110 = axes;
hold on;
plot(T,VnoSin*coV,'k');
%xlabel('Time (s)')
ylabel('Vm (mV)');
title('Membrane potential','fontsize',FS);
set(gca,'ylim',coV*[-0.13 0.04])%
set(gca,'TickDir','out');
set(gca,'fontsize',FS);
set(gca,'ytick',coV*[-0.1 -0.075 -0.05 -0.025 0 0.025])
Cst = stickABC('A',fontA,p1a,p2a,gca);
axb = findobj('Tag','axB')
ax1pos = get(axb,'position');% from the multipanel figure
delete(axb)
newAxcopied = copyobj(ax110,fbig)
set(newAxcopied,'position',ax1pos) %to the
ax100= gca;



figure(11001)
ax11001 = axes;
hold on;
plot(T,Vnoe*coV,'color',[0.5 0.5 0.5]);
xlabel('Time (s)')
ylabel('Vm (mV)');
title('Membrane potential','fontsize',FS);
set(gca,'ylim',coV*[-0.13 0.04])%
set(gca,'TickDir','out');
set(gca,'fontsize',FS);
set(gca,'ytick',coV*[-0.1 -0.075 -0.05 -0.025 0 0.025])
Est = stickABC('E',fontA,p1a,p2a,gca);
axe = findobj('Tag','axNE')
ax1pos = get(axe,'position');% from the multipanel figure
delete(axe)
newAxcopied = copyobj(ax11001,fbig)
set(newAxcopied,'position',ax1pos) %to the
ax1001= gca;



%%% before is when mm == 1 for balance and for no sin;

figure(111)
ax111 = axes;
hold on;
plot(T,V*coV,'k');
%xlabel('Time (s)')
ylabel('Vm (mV)');
title('Membrane potential','fontsize',FS);
set(gca,'ylim',coV*[-0.13 0.04])%
set(gca,'TickDir','out');
set(gca,'fontsize',FS);
set(gca,'ytick',coV*[-0.1 -0.075 -0.05 -0.025 0 0.025])
Cst = stickABC('C',fontA,p1a,p2a,gca);
axc = findobj('Tag','axC')
ax1pos = get(axc,'position');% from the multipanel figure
delete(axc)
newAxcopied = copyobj(ax111,fbig)
set(newAxcopied,'position',ax1pos) %to the
ax1= gca;

%ylim([-1e-12 0]);
%current real
figure(1110)
ax1110 = axes;
plot(T,Iinj*coI,'k');
%plot(T,IIm-(1e-11));
%ylim([-1.5e-9 5e-9])
ylabel('I(nA)')
%xlabel('Time (s)')
title('Current','fontsize',FS);
set(gca,'ylim',coI*[-2e-9 2e-9])
set(gca,'TickDir','out');
set(gca,'fontsize',FS);
Dst = stickABC('D',fontA,p1a,p2a);
box off

axd = findobj('Tag','axD')
ax1pos = get(axd,'position');% from the multipanel figure
delete(axd)
newAxcopied = copyobj(ax1110,fbig)
set(newAxcopied,'position',ax1pos) %to the

%GE and GI
figure(112)
ax112 = axes;
hold on;
plot(T,coG*GE,'g','linewidth',1);
plot(T,coG*GI,'r','linewidth',1);

%plot(T,coG*(GI+GE),'k','linewidth',1);
plot(T,T*0+coG*GL,'k--');
%plot(T,ge,'g')
%plot(T,gi,'r')
xlabel('Time (s)','fontsize',FS);
ylabel('Conductance (nS)','fontsize',FS);
set(gca,'ytick',coG*[0 0.5 1]*1e-7);
set(gca,'TickDir','out');
yl = get(gca,'ylim');
title('Total conductance','fontsize',FS);
legend({'Ge','Gi'})
set(gca,'ytick',[0 10 20 30 40 50]);
set(gca,'ylim',[ -1 50])
set(gca,'fontsize',FS);
Ast = stickABC('B',fontA,p1a,p2a);

axa = findobj('Tag','axA')
ax1pos = get(axa,'position');% from the multipanel figure
delete(axa)
newAxcopied = copyobj(ax112,fbig)
set(newAxcopied,'position',ax1pos) %to the
%filtered voltgages and envelopes. f1

figure(113)
ax113 = axes;
hold on;
plot(T,coV*VF1,'color',colorz1);
hold on;
%plot(T,coV*evf1,'k','linewidth',1);
%xlabel('time (S)','fontsize',FS);
ylabel('Filtered Vm(mV)','fontsize',FS);
%legend({'filtered Vm','Envelope'})
title('BandPass V(f1)','fontsize',FS);
set(gca,'fontsize',FS);
%E1st = stickABC('E1',fontA,p1a,p2a);
E1st = stickABC('F',fontA,p1a,p2a);


axe1 = findobj('Tag','axE1')
ax1pos = get(axe1,'position');% from the multipanel figure
delete(axe1)
newAxcopied = copyobj(ax113,fbig)
set(newAxcopied,'position',ax1pos) %to the
%filtered voltgages and envelopes. f2



figure(114)
ax114 = axes;
hold on;
plot(T,coV*VF2,'color',colorz2);
hold on;
%plot(T,coV*evf2,'k','linewidth',1);
%xlabel('time (s)','fontsize',FS);
ylabel('Filtered Vm(mV)','fontsize',FS);
%legend({'filtered Vm','Envelope'})
title('BandPass V(f2)','fontsize',FS);
%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
set(gca,'fontsize',FS);
%F1st = stickABC('F1',fontA,p1a,p2a);
F1st = stickABC('H',fontA,p1a,p2a);

axf1 = findobj('Tag','axF1')
ax1pos = get(axf1,'position');% from the multipanel figure
delete(axf1)
newAxcopied = copyobj(ax114,fbig)
set(newAxcopied,'position',ax1pos) %to the

figure(115)
ax115 = axes;
hold on;
plot(T,coI*IF1,'color',colorz1);
hold on;
%plot(T,coV*evf2,'k','linewidth',1);
%xlabel('time (s)','fontsize',FS);
ylabel('Filtered I(nA)','fontsize',FS);
%legend({'filtered I','Envelope'})
title('BandPass I(f1)','fontsize',FS);
%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
set(gca,'fontsize',FS);
%E2st = stickABC('E2',fontA,p1a,p2a);
E2st = stickABC('G',fontA,p1a,p2a);

axe2 = findobj('Tag','axE2')
ax1pos = get(axe2,'position');% from the multipanel figure
delete(axe2)
newAxcopied = copyobj(ax115,fbig)
set(newAxcopied,'position',ax1pos) %to the

figure(1155)
ax1155 = axes;
hold on;
plot(T,coI*IF2,'color',colorz2);
hold on;
%plot(T,coV*evf2,'k','linewidth',1);
xlabel('Time (s)','fontsize',FS);
ylabel('Filtered I(nA)','fontsize',FS);
%legend({'filtered I','Envelope'})
title('BandPass I(f2)','fontsize',FS);
%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
set(gca,'fontsize',FS);
%F2st = stickABC('F2',fontA,p1a,p2a);
F2st = stickABC('I',fontA,p1a,p2a);

axf2 = findobj('Tag','axF2')
ax1pos = get(axf2,'position');% from the multipanel figure
delete(axf2)
newAxcopied = copyobj(ax1155,fbig)
set(newAxcopied,'position',ax1pos) %to the

% plot the two Z together and R
figure(116)
ax116 = axes;
hold on;
plot(T,coZ*abs(z1),'linewidth',1,'color',colorz1);
hold on;
plot(T,coZ*abs(z2),'linewidth',1,'color',colorz2);
plot(T,coZ*abs(re),'linewidth',1);
xlabel('Time (s)','fontsize',FS);
ylabel(strcat('Imedance ', '(M','\Omega',')'));
legend({'Imp z(f1)','Imp z(f2)','Rs'},'fontsize',FS);
L = round(length(z1)/4);
mz1 = coZ*mean(real(re(1:L:end-L)));
set(gca,'ylim',[coZ*RR-0.5 coZ*RR+2]);
title('Impedance','fontsize',FS);
set(gca,'fontsize',FS);
%Gst = stickABC('G',fontA,p1a,p2a);
Gst = stickABC('J',fontA,p1a,p2a);

%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
axg = findobj('Tag','axG')
ax1pos = get(axg,'position');% from the multipanel figure
delete(axg)
newAxcopied = copyobj(ax116,fbig)
set(newAxcopied,'position',ax1pos) %to the
legend({'Imp z(f1)','Imp z(f2)','Rs'},'fontsize',5);

% plot the measure cell conductance (gi+ge+gl)


%% impdedance again start figure 3.
fbig2 = open('skel_figure2_method2_nonfilled_f2part2.fig');

% plot the two Z together and R
figure(1166)
ax116 = axes;
hold on;
plot(T,coZ*real(z1),'linewidth',1,'color',colorz1);
hold on;
plot(T,coZ*real(z2),'linewidth',1,'color',colorz2);
plot(T,coZ*real(re),'linewidth',1);
xlabel('Time (s)','fontsize',FS);
ylabel(strcat('Imedance ', '(M','\Omega',')'));
legend({'Imp z(f1)','Imp z(f2)','Rs'},'fontsize',FS);
L = round(length(z1)/4);
mz1 = coZ*mean(real(re(1:L:end-L)));
set(gca,'ylim',[coZ*RR-0.5 coZ*RR+2]);
title('Impedance','fontsize',FS);
set(gca,'fontsize',FS);
%Gst = stickABC('G',fontA,p1a,p2a);
Gst = stickABC('A',fontA,p1a,p2a);

%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
axg2 = findobj('Tag','axG2')
ax1pos = get(axg2,'position');% from the multipanel figure
delete(axg2)
newAxcopied = copyobj(ax116,fbig2)
set(newAxcopied,'position',ax1pos) %to the
legend({'Imp z(f1)','Imp z(f2)','Rs'},'fontsize',5);



figure(120)
ax1120 = axes;
hold on;

hold on;
%plot(T,coG*(GE+GI+GL),'linewidth',1);
plot(T,coG*real(GT),'k','linewidth',0.5);
xlabel('Time (s)','fontsize',FS);
ylabel(strcat('Conductance (nS)'),'fontsize',FS);
legend({'Measured conductance','Real conductance'},'fontsize',FS);
title('Measured Conductane');
stdg = std(coG*real(GT));
set(gca,'ylim',[0 stdg*4]);
set(gca,'fontsize',FS);
%plot(T,(evf1./ampSin1)*MC,'k','linewidth',1);;
Ist = stickABC('C',fontA,p1a,p2a);

axi = findobj('Tag','axI')
ax1pos = get(axi,'position');% from the multipanel figure
delete(axi)
newAxcopied = copyobj(ax1120,fbig2)
set(newAxcopied,'position',ax1pos) %to the
% now plot reconstructed ge gi and

% now plot clean voltage by filtering.(global VC2)
figure(119)
ax119 = axes;
hold on;
plot(T,coV*(VC),'k','linewidth',1);
xlabel('Time (s)','fontsize',FS);
ylabel('Vm(mV)','fontsize',FS);
legend({'Vm, Bandstop for (f1,f2)'},'fontsize',FS);
title('Cleaned Vm','fontsize',FS);
set(gca,'ylim',[-80 -0]);
set(gca,'fontsize',FS);
%Hst = stickABC('H',fontA,p1a,p2a);
Hst = stickABC('B',fontA,p1a,p2a);


axh = findobj('Tag','axH')
ax1pos = get(axh,'position');% from the multipanel figure
delete(axh)
newAxcopied = copyobj(ax119,fbig2)
set(newAxcopied,'position',ax1pos) %to the
% now plot reconstructed ge gi and

figure(121)
ax121 = axes;
hold on;
smpoints = 1;
plot(T,smooth(coG*(real(ge)),smpoints),'g','linewidth',1);
plot(T,smooth(coG*(real(gi)),smpoints),'r','linewidth',1);
plot([T(1) T(end)],coG*[real(gl(1)) real(gl(1))],'k--','linewidth',1);
xlabel('Time (s)','fontsize',FS);
title('Conductanes','fontsize',FS);
stdg = std(coG*real(GT));
set(gca,'ylim',[0 50]);
ylabel(strcat('Conductance (nS)'),'fontsize',FS);
legend({'Ge','Gi'},'fontsize',FS);
title('E and I conductances','fontsize',FS);
set(gca,'ytick',[0 10 20 30 40 50]);
set(gca,'fontsize',FS);
set(gca,'ytick',coG*[0 0.5 1]*1e-7);
set(gca,'ytick',[0 10 20 30 40 50]);
set(gca,'ylim',[ -1 50])
%Jst = stickABC('J',fontA,p1a,p2a);
Jst = stickABC('D',fontA,p1a,p2a);

%plot inset in K


axj = findobj('Tag','axJ')
ax1pos = get(axj,'position');% from the multipanel figure
delete(axj)
newAxcopied = copyobj(ax121,fbig2)
set(newAxcopied,'position',ax1pos) %to the
% now plot reconstructed ge gi and
if 1
    
    figure(1221)
    ax1221 = axes;
    hold on;
    t11 = round(0.9/dt); t22 = round(1.2/dt);
    Tinset = T(t11:t22);
    gein  = real(ge(t11:t22));
    giin  = real(gi(t11:t22));
    GEin  = GE(t11:t22);
    GIin  = GI(t11:t22);
    plot(Tinset,smooth(coG*(real(gein)),smpoints),'g','linewidth',1.5);
    %plot(Tinset,smooth(coG*(real(giin)),smpoints),'r','linewidth',1);
    plot(Tinset,coG*GEin,'color',[0.0 0.7 0.0],'linewidth',0.5);
    %plot(Tinset,coG*GIin,'r--','linewidth',0.5);
    
    box off
    axis off
    set(gca,'ylim',[-10 40]);
    Lst = stickABC('L',fontA,p1a,p2a);
    
    axjins = findobj('Tag','axKinset1')
    ax1pos = get(axjins,'position');% from the multipanel figure
    delete(axjins)
    newAxcopied = copyobj(ax1221,fbig2)
    set(newAxcopied,'position',ax1pos) %to the
    gca;
    lg = legend({'estimated','real'});
    set(lg,'box','off')
    
    
    figure(1222)
    ax1222 = axes;
    hold on;
    t11 = round(0.9/dt); t22 = round(1.2/dt);
    Tinset = T(t11:t22);
    gein  = real(ge(t11:t22));
    giin  = real(gi(t11:t22));
    GEin  = GE(t11:t22);
    GIin  = GI(t11:t22);
    %plot(Tinset,smooth(coG*(real(gein)),smpoints),'g','linewidth',1);
    plot(Tinset,smooth(coG*(real(giin)),smpoints),'r','linewidth',1.5);
    %plot(Tinset,coG*GEin,'g--','linewidth',0.5);
    plot(Tinset,coG*GIin,'color',[0.7 0.0 0.0],'linewidth',0.5);
    pss  =get(gca,'position');
    
    plot([Tinset(end)-0.1 Tinset(end)],[-5  -5],'k') % scale bar of 100 ms.
    plot([Tinset(end) Tinset(end)]+0.01,[-5  5],'k')
    
    box off
    axis off
    set(gca,'position',pss);
    %set(gca,'ylim',ysi1);
    set(gca,'ylim',[-10 40]);
    
    axjins = findobj('Tag','axKinset2')
    ax1pos = get(axjins,'position');% from the multipanel figure
    delete(axjins)
    newAxcopied = copyobj(ax1222,fbig2)
    
    set(newAxcopied,'position',ax1pos) %to the
    gca;
    lg = legend({'estimated','real'});
    set(lg,'box','off')
    plg = get(lg,'position');
    plg(4) = plg(4)*0.6;
    set(lg,'position');
    
    
    
    %%
end


%feb 2023
close all

figure; plot(GE+GI+GL); hold on; plot(ge+gi+gl);
figure; plot(real(z1)); hold on; plot(real(z2))
