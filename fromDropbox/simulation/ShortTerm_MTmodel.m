function [epscamps SPT  RR EE II] = ShortTerm_MTmodel(SPT,u,trec,taui,ase);
% this program will simulate the feedforard thelamic inputs
% SPT are the times of the spikes
%u utilization
%trec tau recovery
%ase determines the size of the first epsc.



epsc(1) = u*ase;
R = 1;
ii = 0;
n = length(SPT);

I = 0;
R = 1;
E = 1-R-I;
k = 0;
epscamps = [];
RR = [];
II = [];
EE = [];
U = u;
epscamps= [];
dt = 0.001; %millisecond by default.
SPTi = round(SPT./dt);
for t = 0:dt:SPT(end)+dt;
    ii = ii+1;
    if ~isempty(find(ii == SPTi))
        PP = 1;
    else
        PP =0;
    end;
    dR = dt*(I/trec)-(PP*R*U); % the second term was added in a correction to the paper in PNAS
    %dR = dt*(I/trec);%-(PP*R*U);
    term2 = PP*U*R;
    dE = dt*(-E/taui)+(term2);
   
    
    E = E+dE;
    R = R+dR;
    I = 1-E-R;
    if PP == 1;
    k = k+1;    
    epscamps(k) = E*ase;
    EE(k) = E;
    II(k)= I;
    RR(k) = R;
    end;
    
end;