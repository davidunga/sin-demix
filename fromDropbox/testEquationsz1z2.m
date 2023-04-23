
%get the measurements from running junk005;
geMeas= ge;
giMeas = gi;
glMeas = gl;
z1Meas = z1;
z2Meas = z2;RE
% now get the measure from the analytical equations without the simulation.
% 
gTT = GL+GI+GE;
w1 = 2*pi*ff;
w2 = 2*pi*ff2;
REE = RE+zeros(1,length(GE));
z1 = REE+1./(gTT+1i*2*pi*ff*c);
p1 = angle(z1);
z2 = REE+1./(gTT+1i*2*pi*ff2*c);
p2 = angle(z2);
z10 = 0+1./(gTT+1i*2*pi*ff*c);
z20 = 0+1./(gTT+1i*2*pi*ff2*c);

ree = REE
c = 1.5e-10;

%test now:
Gtotal = -(1i*(sqrt(c *(z1 - z2).*(w1 - w2).* (c* z1* w1 - c* z2* w1 - c* z1* w2 + c* z2* w2 + 4*1i)) + c*z1* w1 - c* z2* w1 + c* z1* w2 - c* z2* w2))./(2 *(z1 - z2));

Rs_meas = (1/(2*(1i*c*w1 - 1i*c*w2))).*(1i*c*w1*z1 - 1i*c*w2*z1 + 1i*c*w1*z2 - 1i*c*w2*z2 - sqrt((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 +... 
    1i*c*w2*z2).^2 - 4*(1i*c*w1 - 1i*c*w2).*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)));


%Gtotal2 = 1./(z2-Rs_meas)-1i*c*w2;
Gtotal2 = 1./(z1-Rs_meas)-1i*c*w1;
%Gtotal = -(1i*(sqrt(c *(z1 - z2).*(w1 - w2).* (c* z1* w1 - c* z2* w1 - c* z1* w2 + c* z2* w2 + 4*1i)) + c*z1* w1 - c* z2* w1 + c* z1* w2 - c* z2* w2))./(2 *(z1 - z2));
Rs_meas = 1./(-Gtotal2+1i*c*w2)+z2;
Rs_meas2 = Rs_meas;


figure
subplot(2,2,1);
plot(real(Rs_meas)/MC)
hold on;
plot((REE)/MC)
title('Rs real and Rs found, old method');

subplot(2,2,2);
plot(-real(Gtotal))
hold on;
plot(real(gTT)+0.0000000020)
title('G real and G found, old method');
%%%%%%%%%% best from here##############
%Solve[Abs(r + 1/(g + I*w1*c)) == Abs(z1) && Abs(r + 1/(g + I*w2*c)) ==
%Abs(z2), {r, g}] From the paper?
  z1 = (z1); %if I take abs of z1 it does not work. why?
  z2 = (z2);
Re_fromAbs = (1/(2*(1i*c*w1 - 1i*c*w2)))*(1i*c*w1*z1 - 1i*c*w2*z1 + 1i*c*w1*z2 -... 
  1i*c*w2*z2 + ((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 + 1i*c*w2*z2).^2 -... 
     4*(1i*c*w1 - 1i*c*w2)*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)).^0.5);
 
 Gtotal_fromAbs = (1./(-z1 + z2)).*(I*c*w1*z1 + (c^2*w1^2*z1)./(2*(1i*c*w1 - 1i*c*w2)) - (...
    c^2*w1*w2*z1)./(1i*c*w1 - 1i*c*w2) + (c^2*w2^2*z1)./(...
    2*(1i*c*w1 - 1i*c*w2)) - 1i*c*w2*z2 + (c^2*w1^2*z2)./(...
    2*(1i*c*w1 - 1i*c*w2)) - (c^2*w1*w2*z2)./(1i*c*w1 - 1i*c*w2) + (...
    c^2*w2^2*z2)./(2*(1i*c*w1 - 1i*c*w2)) + (1./(2*(1i*c*w1 - 1i*c*w2)))*...
    1i*c*w1* ((1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 +... 
          1i*c*w2*z2).^2 - 4*(1i*c*w1 - 1i*c*w2).*(z1 - z2 + 1i*c*w1*z1.*z2 -... 
           1i*c*w2*z1.* z2)).^0.5 - (1./(2*(1i*c*w1 - 1i*c*w2))).*1i*c*w2*((-1i*c*w1*z1 + 1i*c*w2*z1 - 1i*c*w1*z2 + 1i*c*w2*z2).^2 - ...
        4*(1i*c*w1 - 1i*c*w2).*(z1 - z2 + 1i*c*w1*z1.*z2 - 1i*c*w2*z1.*z2)).^0.5);
 
 

GtotalfromAbs_andimag = 1./(z1-Re_fromAbs)-1i*c*w1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%or:
%Solve[Abs (r + 1/(g + I*w1*c)) == Abs (z1), {g}]
GfromAbs = -((1i*(1i + c*Re_fromAbs*w1 - c*w1*z1))./(Re_fromAbs - z1)); % this is how we use in the findFEGiMultiFrewq_v005_temp function
%%%%%%%%%%%%%%%%%%
%%%%%%%%%% to here best from here##############
r = real(Re_fromAbs);
r = REE;
%r = abs(Re_fromAbs);

GtotalfromAbs_findmethod =(-r - (-(c^2)*(r.^4)*(w1^2) + z1.^2 + (2*(c^2)*(r.^2)*(w1^2)).*z1.^2 - ...
  (c^2)*(w1^2)*(z1.^4)).^0.5)./(r.^2 - z1.^2); 

%did it work?

subplot(2,2,3);
hold on;
plot((REE/MC))
hold on;
plot(real(Re_fromAbs/MC)+0.005,'r'); %add some to see two curves...
set(gca,'ylim',[25 55])
title('R real and R found, new method');
%from theory see 3 line
A = (gTT.^2+(c*w1)^2);
w1 = 2*pi*ff;
ZZ1 = (ree.^2+(2*ree.*gTT)./A+(gTT.^2)./A.^2+((w1*c)^2)./A.^2).^0.5; 

subplot(2,2,4);
hold on;
plot((gTT),'k')
hold on;
%plot(real(GtotalfromAbs_findmethod),'r');
plot(real(GtotalfromAbs_andimag+0.00000002),'g');
plot(real(-GfromAbs+0.000000001),'m'); %good %
%the two give the same!!! 2023
%plot(real(-GfromAbs+0.000000001),'m'); %oct 2022
title('G real and G found, new method, small shift applied to see both');

%%%%%%%% from Gal 2022
  Gtotal= (sqrt(c)*(1i*sqrt(c)*w2*sqrt(w1 - w2).*sqrt(z1 - z2).*sqrt(c*w1.*z1 - c*w1.*z2 - c*w2.*z1 + c*w2.*z2 + 1i*4) + 1i*c*w1*w2.*z1 - 1i*c*w1*w2.*z2 - 1i*c*w2^2.*z1 + 1i*c*w2^2.*z2 + 2*w1 - 2*w2))./(-1*sqrt(w1 - w2).*sqrt(z1 - z2).*sqrt(c*w1.*z1 - c*w1.*z2 - c*w2.*z1 + c*w2.*z2 + 1i*4) - sqrt(c)*w1.*z1 + sqrt(c)*w1.*z2 + sqrt(c)*w2.*z1 - sqrt(c)*w2.*z2);
  
  
