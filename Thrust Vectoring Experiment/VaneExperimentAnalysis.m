%calibration
clc
clear all
close all
raw=readmatrix("GMotorTest.xlsx")
time=raw(:,1);
trval=raw(:,2);
tqP=raw(:,3);
tq1P=raw(:,4);
tq2P=raw(:,5);
normalizedthrust=trval*17/max(trval)


% figure
% plot(time,trval)
% xlim([0.615752,0.615779])
% ylim([0,20])
% grid on
% xlabel('UTC time')
% ylabel('Thrust (whatever the output is)')
% title('G moto booyah')


% [inpict map] = imread('thrustGcurve.png');
% inpict = ind2rgb(inpict,map); % this is an indexed image.  it must be RGB
% inpict = flipud(inpict);
% hold on
% plot(time,normalizedthrust)
% xlim([0.615752,0.615779])
% grid on
% xlabel('UTC time')
% ylabel('Thrust (N)')
% title('G motor Thust Curve')
% hi = image(inpict,'xdata',[-1 1],'ydata',[0 3]);
% set(gca,'ydir','normal')
% uistack(hi,'down') % move it beneath the scatter plot

% Load the image
[inpict, map] = imread('Thrug.png');

% Convert indexed image to RGB if it's indexed
if ~isempty(map)
    inpict = ind2rgb(inpict, map);
end

% Flip the image vertically if needed
inpict = flipud(inpict);

% Create a figure and hold for multiple plots
figure;
hold on;

% Get the current axis limits for proper image placement
xLimits = [0.615753, 0.615778];
yLimits = [0, 20];

% Place the image in the background with proper alignment to axis limits
hi = image(xLimits, yLimits, inpict);
uistack(hi, 'bottom');  % Move image to bottom layer

% Plot your data
plot(time, normalizedthrust,'k-','LineWidth',2);
plot([1 1],[1 1],'r-','LineWidth',2)
xlim(xLimits);
ylim(yLimits);  % Ensure y-axis matches the image bounds

% Add grid and labels
grid on;
xlabel('UTC time');
ylabel('Thrust (N)');
title('Thrust Curve Calibration');
legend('Experimental Data','G75J-10A Reference Curve')
% Ensure correct y-direction
set(gca, 'YDir', 'normal');
%%
figure;
hold on;
plot(time, tqP, 'LineWidth', 1.5);
plot(time, tq1P, 'LineWidth', 1.5);
plot(time, tq2P, 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Torque');
title('Torque Comparison');
legend('tqP', 'tqP1', 'tqP2');
hold off;
%% H motor run1 
clear
raw1=readmatrix("HMotorTest1.csv");
raw2=readmatrix("HMotorTest2Actual.csv");
raw3=readmatrix("HMotorTest3.csv");
raw4=readmatrix("HMotorTest4Actual.csv");
raw1=raw1(742:948,:);
raw2=raw2(111:315,:);
raw3=raw3(131:335,:);
raw4=raw4(192:400,:);
% time=[0:0.1:20];
n = length(raw1);

time = linspace(0, (n-1)*0.1, n)';
trval=raw1(:,2); tqP1=raw1(:,3); tq1P1=raw1(:,4); tq2P1=raw1(:,5);
trval2=raw2(:,2); tqP2=raw2(:,3); tq1P2=raw2(:,4); tq2P2=raw2(:,5);
trval3=raw3(:,2); tqP3=raw3(:,3); tq1P3=raw3(:,4); tq2P3=raw3(:,5);
trval4=raw4(:,2); tqP4=raw4(:,3); tq1P4=raw4(:,4); tq2P4=raw4(:,5);
% %1
figure
plot(time,trval,"LineWidth",2)
xlabel('Time (S)')
ylabel('Thrust (Loadcell Reading)')
title('Thrust vs Time Run 1')
xlim([0,20])
ylim([0,4000])
figure
plot(time,tqP1,time,tq1P1,time,tq2P1,"LineWidth",2)
title('Moment vs time Run1')
xlabel('Time(S)')
ylabel('Moment (N-M)')
xlim([0,20])

%2
n = length(raw2);
time2 = linspace(0, (n-1)*0.1, n)';

figure
plot(time2,trval2,"LineWidth",2)
xlabel('Time (S)')
ylabel('Thrust (Loadcell Reading)')
title('Thrust vs Time Run 2')
xlim([0,20])
ylim([0,4000])
figure
plot(time2,tqP2,time2,tq1P2,time2,tq2P2,"LineWidth",2)
title('Moment vs time Run 2')
xlabel('Time(S)')
ylabel('Moment (N-M)')
xlim([0,20])

%3
n = length(raw3);
time3 = linspace(0, (n-1)*0.1, n)';


figure
plot(time3,trval3,"LineWidth",2)
xlabel('Time (S)')
ylabel('Thrust (Loadcell Reading)')
title('Thrust vs Time Run 3')
xlim([0,20])
ylim([0,4000])
figure
plot(time3,tqP3,time3,tq1P3,time3,tq2P3,"LineWidth",2)
title('Moment vs time Run 3')
xlabel('Time (S)')
ylabel('Moment (N-M)')
xlim([0,20])

%4
n = length(raw4);
time4 = linspace(0, (n-1)*0.1, n)';


figure
plot(time4,trval4,"LineWidth",2)
xlabel('Time (S)')
ylabel('Thrust (Loadcell Reading)')
title('Thrust vs Time Run 4')
xlim([0,20])
ylim([0,4000])
figure
plot(time4,tqP4,time4,tq1P4,time4,tq2P4,"LineWidth",2)
title('Moment vs time Run 4')
xlabel('Time (S)')
ylabel('Moment (N-M)')
xlim([0,20])

%%
%all thrust curves
figure
hold on
rgb1=[101, 111, 122]/255
rgb2= [167, 177, 185]/255
rgb3= [176, 0, 0]/255
rgb4=[240, 64, 69]/255
plot(time, trval, "LineWidth", 1.5, "Color", rgb1)     % Orange
plot(time2, trval2, "LineWidth", 1.5, "Color", rgb2)        % Blue
plot(time3, trval3, "LineWidth", 1.5, "Color", rgb3)   % Green
plot(time4, trval4, "LineWidth", 1.5, "Color", rgb4)   % Dark Red
grid on
xlabel('$\mathrm{Time}\ (s)$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\mathrm{Force}\ (\mathrm{Load\ Cell\ Reading})$', 'Interpreter', 'latex', 'FontSize', 12)
title('$\mathrm{Force\ vs\ Time\ H\ motors}$', 'Interpreter', 'latex', 'FontSize', 14)
legend('Run 1', 'Run 2', 'Run 3', 'Run 4')
set(gca, 'TickLabelInterpreter', 'latex')
xlim([0,20])
ylim([0,4000])
%all moments
figure
hold on

% Run 1 - Orange color group
h1 = plot(time, tqP1, "LineWidth", 1.5, "Color", [111, 167, 237]/255);
h2 = plot(time, tq1P1, "LineWidth", 1.5, "Color", [79, 104, 137, 0.7*255]/255, "LineStyle", ":");
h3 = plot(time, tq2P1, "LineWidth", 1.5, "Color", [121, 148, 186, 0.7*255]/255, "LineStyle", ":");

% Run 2 - Blue color group
h4 = plot(time2, tqP2, "LineWidth", 1.5, "Color", [125, 131, 139]/255);
h5 = plot(time2, tq1P2, "LineWidth", 1.5, "Color", [83, 91, 99, 0.7*255]/255, "LineStyle", ":");
h6 = plot(time2, tq2P2, "LineWidth", 1.5, "Color", [195, 199, 201, 0.7*255]/255, "LineStyle", ":");

% Run 3 - Green color group
h7 = plot(time3, tqP3, "LineWidth", 1.5, "Color", [254, 75, 80]/255);
h8 = plot(time3, tq1P3, "LineWidth", 1.5, "Color", [144, 105, 106, 0.7*255]/255, "LineStyle", ":");
h9 = plot(time3, tq2P3, "LineWidth", 1.5, "Color", [191, 148, 149, 0.7*255]/255, "LineStyle", ":");

% Run 4 - Dark Red color group
h10 = plot(time4, tqP4, "LineWidth", 1.5, "Color",[247, 228, 95]/255);
h11 = plot(time4, tq1P4, "LineWidth", 1.5, "Color", [150, 122, 103, 0.7*255]/255, "LineStyle", ":");
h12 = plot(time4, tq2P4, "LineWidth", 1.5, "Color", [200, 164, 133, 0.7*255]/255, "LineStyle", ":");

% LaTeX formatting for labels
title('$\mathrm{Moment\ vs\ Time}$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\mathrm{Time}\ (s)$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\mathrm{Moment}\ (N\mbox{-}m)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex')

% Add legend with LaTeX formatting - rearranged to have all Average lines first
legend([h1, h4, h7, h10, h2, h5, h8, h11, h3, h6, h9, h12], ...
       {'Run 1 - Average', 'Run 2 - Average', 'Run 3 - Average', 'Run 4 - Average', ...
        'Run 1 - LC 1', 'Run 2 - LC 1', 'Run 3 - LC 1', 'Run 4 - LC 1', ...
        'Run 1 - LC 2', 'Run 2 - LC 2', 'Run 3 - LC 2', 'Run 4 - LC 2'}, ...
       'Interpreter', 'latex', 'Location', 'best')
       
xlim([0,20])

%% moment avg main
figure
hold on

% Run 1 - Orange color group
h1 = plot(time, tqP1, "LineWidth", 1.5, "Color", [111, 167, 237]/255);

% Run 2 - Blue color group
h4 = plot(time2, tqP2, "LineWidth", 1.5, "Color", [125, 131, 139]/255);


% Run 3 - Green color group
h7 = plot(time3, tqP3, "LineWidth", 1.5, "Color", [254, 75, 80]/255);


% Run 4 - Dark Red color group
h10 = plot(time4, tqP4, "LineWidth", 1.5, "Color",[247, 228, 95]/255);


% LaTeX formatting for labels
title('$\mathrm{Moment\ vs\ Time}$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\mathrm{Time}\ (s)$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\mathrm{Moment}\ (N\mbox{-}m)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex')
legend('Run 1 - Average', 'Run 2 - Average', 'Run 3 - Average', 'Run 4 - Average','Interpreter', 'latex', 'Location', 'best')
xlim([0,20])

%% Servo
timeservo = [0, 1:1:4, 1004:1:1007, 2007:1:2010, 3010:1:3013, 4013:1:4016, 5016, 5017:1:5031, ...
        6031, 6032:1:6034, 7034, 7035:1:7037, 8037, 8038:1:8040, 9040, 9041:1:9043, ...
        10043, 10044:1:10046, 11046, 11047:1:11061, 21061];

% Convert time to seconds for plotting
time_sec = timeservo/1000;

% Servo angle data
servo_angle = [0, 90:1:93, 93, 94:1:96, 96, 97:1:99, 99, 100:1:102, 102, 103:1:105, ...
              105, 104:-1:90, 90, 89:-1:87, 87, 86:-1:84, 84, 83:-1:81, 81, 80:-1:78, ...
              78, 77:-1:75, 75, 76:1:90, 90];

% Calculate AoA
aoa = servo_angle - 90;
aoa(1)=0;

% Plot AoA vs Time


figure
plot(time_sec, aoa, 'LineWidth', 1.5)
grid on
title('$\mathrm{Angle\ of\ Attack\ vs\ Time}$', 'Interpreter', 'latex', 'FontSize', 14)
xlabel('$\mathrm{Time}\ (s)$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\mathrm{Angle\ of\ Attack}\ (^\circ$)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex')
xlim([0, 22])
%% moment run 3 with servo. 

%%%%%%%%%%%%%% imperial conversion
% Scale the torque values
tqP3 = tqP3 * 0.7376;
tq1P3 = tq1P3 * 0.7376;
tq2P3 = tq2P3 * 0.7376;
%%
% Create a tiled layout
t=tiledlayout(1, 2);

% Set default text interpreter to LaTeX for all text
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesFontName', 'LaTeX');
set(0, 'DefaultTextFontName', 'LaTeX');

% First tile - Original plot
nexttile;
plot(time_sec, aoa, 'LineWidth', 1.5, "Color", [247, 228, 95]/255)
hold on
plot(time3, tqP3, "LineWidth", 1.5, "Color", [167, 177, 185]/255)
plot(time3, tq1P3, "LineWidth", 1.5, "Color", [153, 0, 0]/255)
plot(time3, tq2P3, "LineWidth", 1.5, "Color", [240, 64, 69]/255)
ylim([-20, 20])
% Normalized to max value from thrustcurve.com slightly different burn time
plot(time3, (trval3*45/max(trval))*0.224808, "LineWidth", 1.5, "Color", [111, 167, 237]/255)

% Add vertical lines for steady state operation region
yLimits = ylim;
plot([7 7], yLimits, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot([11 11], yLimits, 'k--', 'LineWidth', 1.5); 

% Add labels and legend
title('\textbf{Run 3: All Measurements}')
xlabel('\textbf{Time (s)}')
ylabel('\textbf{lb-ft $\quad | \quad$ lb $\quad | \quad$ Degrees}')
legend({'Vane Angle of Attack', 'Average Moment', 'Moment Load Cell 1', 'Moment Load Cell 2', 'Axial Thrust', 'Steady State Operation Region'})
grid on  
xlim([0, 20])

% Second tile - AOA vs Moment with linear fit
nexttile;

% Prepare data
tr3M = tqP3(71:110);
tr3M1 = tq1P3(71:110);
tr3M2 = tq2P3(71:110);
trAOA = abs(aoa(42:58));

% Interpolate data
t1 = linspace(0, 1, length(trAOA));
t2 = linspace(0, 1, length(tr3M));
AOA_interp = interp1(t1, trAOA, t2);

% Plot data points
plot(AOA_interp, tr3M2, 'r>', 'DisplayName', 'Load Cell 2')
hold on

% Calculate and add linear fit
a = AOA_interp(:) \ tr3M2(:);  % Least-squares fit with no intercept
x_fit = linspace(0, max(AOA_interp), 100);
y_fit = a * x_fit;
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Fit: y = %.4fx', a))

% Add labels and legend
xlim([0, 15])
ylim([0,3])
xlabel('\textbf{Angle of Attack (Degrees)}', 'Interpreter', 'latex')
ylabel('\textbf{Moment (lb-ft)}', 'Interpreter', 'latex')
title('\textbf{Load Cell 2 vs. Angle of Attack}', 'Interpreter', 'latex')
legend('Location', 'northwest')
grid on

% Adjust overall layout
title(t, '\textbf{Jet Vane Thrust Vector Control Experimental Data}', 'Interpreter', 'latex')
