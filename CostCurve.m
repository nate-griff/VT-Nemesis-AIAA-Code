

%Cost Curve Estimation Via Tactical Missile Design 
%Cost of First Interceptor
C1=12800;
%manufacturing cost at scale
j=900;

%motor cost at scale
m=5500;
C1=C1-j-m;
%Cost Curve Coeficcient
L=0.85;

x=[1:1:10050];

%m=5500*L

cx=C1*L.^(log2(x));
%Adding motor, manufacturing, and recurring costs post curve adjustment
cx=cx+j+m+1400;

%Mean Interceptor Cost only
PlainMean=mean(cx)

%adding 10% Profit Margin
cxx=cx;
margin=cx*0.1;
cx=cx+margin;

fprintf('Mean Missie Cost')
Cmean=mean(cx)

fprintf('Total Manufacturing Cost')
tcost=sum(cx)

fprintf('Missile Cost with Recuring Costs')
%mcrc=Cmean+1400


fprintf('Total Program Cost')
progcost=tcost+1050*(1400)


%fprintf('Cost Percentages')
%PercCostAdmin=(700/mcrc)*100
%PercCostIntegration=(340/mcrc)*100
%RocketPer=(Cmean/mcrc)*100

%Sum of total profit margin
margint=sum(margin)

close all
figure
plot(cxx,LineWidth=1.5)
fontsize(14, 'points')
title('Nemesis Production Cost Curve','FontSize',20, "LineWidth",2)
ylabel('Cost USD FY2024','FontSize',18)
xlabel('Produciton Number','FontSize',18)
set(gca, 'FontName', 'Acherus Grotesque Regular')
ax = gca;
ax.YRuler.Exponent = 0;
grid on 


