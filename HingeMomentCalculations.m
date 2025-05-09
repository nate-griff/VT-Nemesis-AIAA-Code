%Hinge Moment Calculations
%Coefficient Values taken from the report Prediciton of Hinge Moment
%Coefficients for Nose Mounted Canard Controls at Supersonic Seeds by
%Michael G. Landers

%air speed
v2=343*2;
v125=343*1.25;
v15=343*1.5;
%Density at sealevel: Where it is most dense
rho=1.225;
%Control surface area in m2
S=0.00827;
c=.03;


%Hinge Moment Coefficients by aoa
HCM=[0.0, 0;
     0.5, -0.00001467;
     1.0, -0.00002939;
     1.5, -0.00004409;
     2.0, -0.00005879;
     2.5, -0.00007348;
     3.0, -0.00008818;
     3.5, -0.00010287;
     4.0, -0.00011757;
     4.5, -0.00013227;
     5.0, -0.00014696;
     5.5, -0.00016166;
     6.0, -0.00017636;
     6.5, -0.00019105;
     7.0, -0.00020575;
     7.5, -0.00022045;
     8.0, -0.00023514;
     8.5, -0.00024984;
     9.0, -0.00026453;
     9.5, -0.00027923;
    10.0, -0.00029393;
    10.5, -0.00030862];

%Hinge Moment at Mach 1.25
HCM125=[0.0, 0;
     0.5, -0.00002817;
     1.0, -0.00005172;
     1.5, -0.00007758;
     2.0, -0.00010344;
     2.5, -0.00012930;
     3.0, -0.00015516;
     3.5, -0.00018102;
     4.0, -0.00020688;
     4.5, -0.00023274;
     5.0, -0.00025860;
     5.5, -0.00028446;
     6.0, -0.00031032;
     6.5, -0.00033618;
     7.0, -0.00036204;
     7.5, -0.00038790;
     8.0, -0.00041376;
     8.5, -0.00043962;
     9.0, -0.00046548;
     9.5, -0.00049134;
    10.0, -0.00051720;
    10.5, -0.00054306];
%Hinge Moment Coeficient at mach 1.5
HCM15=[0.0, 0;
     0.5, -0.00002935;
     1.0, -0.00005856;
     1.5, -0.00008784;
     2.0, -0.00011711;
     2.5, -0.00014639;
     3.0, -0.00017567;
     3.5, -0.00020495;
     4.0, -0.00023423;
     4.5, -0.00026351;
     5.0, -0.00029279;
     5.5, -0.00032207;
     6.0, -0.00035134;
     6.5, -0.00038062;
     7.0, -0.00040990;
     7.5, -0.00043918;
     8.0, -0.00046846;
     8.5, -0.00049774;
     9.0, -0.00052702;
     9.5, -0.00055629;
    10.0, -0.00058557;
    10.5, -0.00061485];

%Hinge Moment 
HM=0.5*rho*(v2^2)*HCM(:,2)*S*c;
HM125=0.5*rho*(v125^2)*HCM125(:,2)*S*c;
HM15=0.5*rho*(v15^2)*HCM15(:,2)*S*c;
%Hinge Moment Pound 
HMPI=-1*(HM*8.85);
HMPI125=-1*(HM125*8.85);
HMPI15=-1*(HM15*8.85);
close all
%Maybe solve for normal force next.
figure('name','Hinge Moment By Control Surface AOA at Sea Level')
hold on
plot(HCM(:,1),HMPI125)
plot(HCM(:,1),HMPI15)
plot(HCM(:,1),HMPI)
title('Hinge Moment by AOA at Sea Level')
ylabel('Hinge Moment in Inch Pounds')
xlabel('AOA of Control Surface in Degrees')
legend({'Mach 1.25', 'Mach 1.5','Mach 2'},'Location','southeast')
hold off



%Normal Force (Lift) Calculations for the control surfaces.
%Max Turn rate will be from two fins moving at the same time.
%Mach 1.25 Normal Force Coefficient
CNF125=[0.0, 0;
     0.5, 0.00361731;
     1.0, 0.00722098;
     1.5, 0.01083147;
     2.0, 0.01444195;
     2.5, 0.01805244;
     3.0, 0.02166293;
     3.5, 0.02527342;
     4.0, 0.02888391;
     4.5, 0.03249440;
     5.0, 0.03610488;
     5.5, 0.03971537;
     6.0, 0.04332586;
     6.5, 0.04693635;
     7.0, 0.05054684;
     7.5, 0.05415732;
     8.0, 0.05776781;
     8.5, 0.06137830;
     9.0, 0.06498879;
     9.5, 0.06859928;
    10.0, 0.07220977;
    10.5, 0.07582025];
%Norrmal Force at Mach 1.5
CNF15=[0.0, 0;
     0.5, 0.00288718;
     1.0, 0.00577648;
     1.5, 0.00866472;
     2.0, 0.01155297;
     2.5, 0.01444121;
     3.0, 0.01732945;
     3.5, 0.02021769;
     4.0, 0.02310593;
     4.5, 0.02599417;
     5.0, 0.02888241;
     5.5, 0.03177065;
     6.0, 0.03465889;
     6.5, 0.03754714;
     7.0, 0.04043538;
     7.5, 0.04332362;
     8.0, 0.04621186;
     8.5, 0.04910010;
     9.0, 0.05198834;
     9.5, 0.05487658;
    10.0, 0.05776483;
    10.5, 0.06065307];
%Normal Force at Mach 2
CNF2=[0.0, 0;
     0.5, 0.00204116;
     1.0, 0.00408323;
     1.5, 0.00612484;
     2.0, 0.00816646;
     2.5, 0.01020807;
     3.0, 0.01224969;
     3.5, 0.01429130;
     4.0, 0.01633292;
     4.5, 0.01837453;
     5.0, 0.02041615;
     5.5, 0.02245776;
     6.0, 0.02449938;
     6.5, 0.02654099;
     7.0, 0.02858260;
     7.5, 0.03062422;
     8.0, 0.03266583;
     8.5, 0.03470745;
     9.0, 0.03674906;
     9.5, 0.03879068;
    10.0, 0.04083229;
    10.5, 0.04287391];

L125=0.5*rho*(v125)^2*S.*CNF125(:,2);
L15=0.5*rho*(v15)^2*S.*CNF15(:,2);
L2=0.5*rho*(v2)^2*S.*CNF2(:,2);

figure('name','Lift by AOA in Newtons')
hold on 
plot(HCM(:,1),L125)
plot(HCM(:,1),L15)
plot(HCM(:,1),L2)
title('Single Control Surface Lift by AOA')
ylabel('Lift in Newtons')
xlabel('AOA of Control Surface in Degrees')
legend({'Mach 1.25', 'Mach 1.5','Mach 2'},'Location','southeast')
hold off

L125P=0.224809*L125;
L15P=0.224809*L15;
L2P=0.224809*L2;
figure('name','Lift by AOA in Pounds Force')
hold on
plot(HCM(:,1),L125P)
plot(HCM(:,1),L15P)
plot(HCM(:,1),L2P)
title('Single Control Surface Lift by AOA')
ylabel('Lift in Pounds Force')
xlabel('AOA of Control Surface in Degrees')
legend({'Mach 1.25', 'Mach 1.5','Mach 2'},'Location','southeast')
hold off

%%Max Turn Rate Calculation Section
%Cd=0.4;
%predicted Load Factor
%all values in metric

%Extrapolating to 22 deg aoa for lift
C125poly=polyfit(CNF125(:,1),CNF125(:,2),1)
C15poly=polyfit(CNF125(:,1),CNF15(:,2),1);
C2poly=polyfit(CNF125(:,1),CNF2(:,2),1);

maxLM125=C125poly(1)*22
maxLM15=C15poly(1)*22
maxLM2=C2poly(1)*22

%Extrapolating hinge force at 22 deg
M125poly=polyfit(CNF125(:,1),HMPI125,1);
M15poly=polyfit(CNF125(:,1),HMPI15,1);
M2poly=polyfit(CNF125(:,1),HMPI,1);
maxHM125=M125poly(1)*22
maxHM15=M15poly(1)*22
maxHM2=M2poly(1)*22


%Turn Rate
n=20;
Tav=7000; %*0.224809
%mach 2 cl at 10 deg
clmax=0.0429;
S2=2*S;
radsec=((rho*S2*clmax)/(2*65*9.81))^(1/2)*((n^2-1)/n)^(1/2)
rmin=(2*65*9.81)/(rho*9.81*S2*clmax)*(n/(n^2-1)^1/2)
%0.21127706036962632 deg sec



