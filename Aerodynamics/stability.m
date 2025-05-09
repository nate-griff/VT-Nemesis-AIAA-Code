% filepath: /Users/nathangriffin/Downloads/mlab-missile/stability.m

% Load the data from stabmarg.csv
data = readmatrix('stab-final.csv');
time = data(:, 1); % Time (s)
stab_margin = data(:, 2); % Stability Margin

% Filter data to include only time <= 45.045 seconds
valid_indices = time <= 44.045;
time = time(valid_indices);
stab_margin = stab_margin(valid_indices);

% Set default font properties
set(groot, 'DefaultAxesFontName', 'Acherus Grotesque');
set(groot, 'DefaultTextFontName', 'Acherus Grotesque');
set(groot, 'DefaultLegendFontName', 'Acherus Grotesque');

% Create the plot
figure;
plot(time, stab_margin, 'LineWidth', 2.5, 'Color', [142, 17, 23]/255); % Blue color

% Add grid, labels, and title
grid on;
grid minor;
xlabel('Time (s)', 'Interpreter', 'tex', 'FontSize', 20);
ylabel('Stability Margin', 'Interpreter', 'tex', 'FontSize', 20);
% title('Stability Margin Until Apogee', 'Interpreter', 'tex', 'FontSize', 30);

% Adjust axis properties
set(gca, 'FontSize', 18, 'FontName', 'Acherus Grotesque');
xlim([min(time), max(time)]);
ylim([min(stab_margin), max(stab_margin)]);

% Export the plot as a high-resolution image
exportgraphics(gcf, 'stability_margin_plot.png', 'Resolution', 600);