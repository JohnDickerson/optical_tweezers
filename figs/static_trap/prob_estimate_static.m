function [] = prob_estimate_static(fname)
% file containing [y-distance-to-laser z-distance-to-laser prob-of-trap]
%path = 'run_01.txt';
path = fname;%/home/rob/SoftwareStaging/OTRegion1/build/new_probs.txt'; %tweezers/build/bin/new_output.txt';%trapping_probabilities.txt';

% loads in lines of three floats, ignoring lines beginning with #
[y_dist z_dist p_trap] = textread(path, '%f %f %f', 'commentstyle', 'shell');

% sets the grid size for the contour plot
y_vals = [ min(y_dist):2.5e-7:max(y_dist) ] * 1e6;
z_vals = [ min(z_dist):2.5e-7:max(z_dist) ] * 1e6;
min(y_dist)
max(y_dist)
min(z_dist)
max(z_dist)
%y_vals = [ min(y_dist):1e-6:max(y_dist) ] * 1e6;
%z_vals = [ min(z_dist):1e-6:max(z_dist) ] * 1e6;
% force contours at 0.0, 0.1, ..., 0.9+
line_values = 0:0.1:0.9;   
length(p_trap)
length(z_vals)
length(y_vals)
[B,IX] = sort(y_dist,1);
p_trap = p_trap(IX);
% reshape single column of trapping probabilities into a grid
c_matrix = reshape(p_trap, length(z_vals), length(y_vals));

fontsize = 20;
% plots the contours prescribed by line_values, flipping the Y-axis
xlabel('Y (in \mum)','FontSize',fontsize);
ylabel('Z (in \mum)','FontSize',fontsize);
title('Trapping probability for static trap','FontSize',fontsize);
axis([min(y_vals) max(y_vals) -9.0 max(z_vals)])
%title('Trapping probability for moving trap (0.65\mum/s in Y)','FontSize',fontsize);%, 'Interpreter', 'LaTex')
hold on;
contourf(y_vals, z_vals, c_matrix, line_values,'LineStyle',':','LineWidth',1e-100);
set(gca, 'FontSize',fontsize);
set(gca, 'YDir', 'reverse')
hold off;
