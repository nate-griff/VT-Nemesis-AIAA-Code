g=9.81;
mg=3;
close all
v=(0:34.3:1029)
turnradius=(v.^2)./(mg*g)
vv=[-flip(v) v]
turnrad=[flip(turnradius) turnradius]
figure
hold on
title('Target Max G Turn Radius')
plot(v,turnradius)
ylabel('Turn Radius (m)')
xlabel('Velocity m/s')
hold off

maxalt=9144;

%title('Manuvering Range')
%plot(vv,turnrad)
%camroll(-90)
%set(gca, 'ydir', 'reverse')

vlock=linspace(167, 1029, 6);
trl=(vlock.^2)./(mg*g)
t=linspace(0, pi/2, 1000);

[x,y]=pol2cart(t, trl(1));
[x2,y2]=pol2cart(t, trl(2));
[x3,y3]=pol2cart(t, trl(3));
[x4,y4]=pol2cart(t, trl(4));
[x5,y5]=pol2cart(t, trl(5));
[x6,y6]=pol2cart(t, trl(6));
mkm=8.04672;
ftkm=9.144;
figure
hold on
plot(x/1000,y/1000)
plot(x2/1000,y2/1000)
plot(x3/1000,y3/1000)
plot(x4/1000,y4/1000)
plot(x5/1000,y5/1000)
plot(x6/1000,y6/1000)
title('Turn Radius by Mach Number')
xlabel('Turn Radius (km)')
ylabel('Turn Radius (km)')
legend('Mach 0.5' ,'Mach 1', 'Mach 1.5', 'Mach 2', 'Mach 2.5', 'Mach 3')

hold off
%radius more zoomed in on things within the kill area
%irrelevent if it can turn 90 in the box
%ground covered by slight turns within box
%angle covered by object turning burn 10 sec
%or by the time it leaves our area
%%


figure
hold on 
plot(x/1000,y/1000)
plot(x2/1000,y2/1000)
plot(x3/1000,y3/1000)
plot(x4/1000,y4/1000)
plot(x5/1000,y5/1000)
plot(x6/1000,y6/1000)
rectangle('Position',[0 0 mkm mkm])
text(2,5,'5 Mile')
rectangle('Position',[0 0 ftkm ftkm])
text(ftkm+0.5,8,'30,000 ft')
title('Turn Radius by Mach Number')
xlabel('Turn Radius (km)')
ylabel('Turn Radius (km)')
legend('Mach 0.5' ,'Mach 1', 'Mach 1.5', 'Mach 2', 'Mach 2.5', 'Mach 3')

hold off


%turn radius which completes 90deg turn within operational area
voparea=linspace(10, 514.5, 5);
troparea=(voparea.^2)./(mg*g)

[x,y]=pol2cart(t, troparea(1));
[x2,y2]=pol2cart(t, troparea(2));
[x3,y3]=pol2cart(t, troparea(3));
[x4,y4]=pol2cart(t,troparea(4));
[x5,y5]=pol2cart(t, troparea(5));


figure 
hold on 
plot(x/1000,y/1000)
plot(x2/1000,y2/1000)
plot(x3/1000,y3/1000)
plot(x4/1000,y4/1000)
plot(x5/1000,y5/1000)
rectangle('Position',[0 0 mkm mkm])
text(2,5,'5 Mile')
rectangle('Position',[0 0 ftkm ftkm])
text(ftkm+0.5,8,'30,000 ft')
title('Turn Radius within Operatinal Area')
xlabel('Turn Radius (km)')
ylabel('Turn Radius (km)')
legend('10 ms', '136 ms', '262 ms' ,'388 ms', '514 ms')
hold off

[x,y]=pol2cart(t, trl(1));
[x2,y2]=pol2cart(t, trl(2));
[x3,y3]=pol2cart(t, trl(3));
[x4,y4]=pol2cart(t, trl(4));
[x5,y5]=pol2cart(t, trl(5));
[x6,y6]=pol2cart(t, trl(6));
