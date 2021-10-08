%function [] = prob_estimate_static(fname)
% file containing [y-distance-to-laser z-distance-to-laser prob-of-trap]
%path = 'run_01.txt';
%path = fname;%/home/rob/SoftwareStaging/OTRegion1/build/new_probs.txt'; %tweezers/build/bin/new_output.txt';%trapping_probabilities.txt';
clear all;
close all;
clc;

%pathCPU = 'outputCPUFull.txt';
%pathGPU = 'outputGPUFull.txt';
pathCPU = 'final1024CPU.txt';
pathGPU = 'final1024GPU.txt';
%pathCPU = 'Results/outputGPU100.txt';
%pathGPU = 'Results/final1024GPU.txt';



% loads in lines of three floats, ignoring lines beginning with #
[y_distCPU z_distCPU p_trapCPU] = textread(pathCPU, '%f %f %f', 'commentstyle', 'shell');

% sets the grid size for the contour plot
% y_valsCPU = [ min(y_distCPU):2.5e-7:max(y_distCPU) ] * 1e6;
% z_valsCPU = [ min(z_distCPU):2.5e-7:max(z_distCPU) ] * 1e6;
[min(y_distCPU) max(y_distCPU) min(z_distCPU) max(z_distCPU)]

y_valsCPU = [ 0:2.5e-7:20.99e-6] * 1e6;
z_valsCPU = [-20.0e-6:2.5e-7:8.1e-6] * 1e6;

%y_valsCPU = [ min(y_distCPU):1e-6:max(y_distCPU) ] * 1e6;
%z_valsCPU = [ min(z_distCPU):1e-6:max(z_distCPU) ] * 1e6;
[min(y_distCPU) max(y_distCPU) min(z_distCPU) max(z_distCPU)]
% force contours at 0.0, 0.1, ..., 0.9+
line_valuesCPU = 0:0.1:0.9;   
[length(p_trapCPU) length(z_valsCPU) length(y_valsCPU) (length(z_valsCPU)*length(y_valsCPU)) length(p_trapCPU)-(length(z_valsCPU)*length(y_valsCPU))]
[B,IX] = sort(y_distCPU,1);
p_trapCPU = p_trapCPU(IX);
% reshape single column of trapping probabilities into a grid
c_matrixCPU = reshape(p_trapCPU(1:length(z_valsCPU)*length(y_valsCPU)), length(z_valsCPU), length(y_valsCPU));

figure,
fontsize = 20;
% plots the contours prescribed by line_values, flipping the Y-axis
xlabel('Y (in \mum)','FontSize',fontsize);
ylabel('Z (in \mum)','FontSize',fontsize);
title('Trapping probability for static trap','FontSize',fontsize);
axis([min(y_valsCPU) max(y_valsCPU) -9.0 max(z_valsCPU)])
%title('Trapping probability for moving trap (0.65\mum/s in Y)','FontSize',fontsize);%, 'Interpreter', 'LaTex')
hold on;
contourf(y_valsCPU, z_valsCPU, c_matrixCPU, line_valuesCPU,'LineStyle',':','LineWidth',1e-100);
set(gca, 'FontSize',fontsize);
set(gca, 'YDir', 'reverse')
hold off;


% loads in lines of three floats, ignoring lines beginning with #
[y_distGPU z_distGPU p_trapGPU] = textread(pathGPU, '%f %f %f', 'commentstyle', 'shell');

% sets the grid size for the contour plot
%y_valsGPU = [ min(y_distGPU):2.5e-7:max(y_distGPU) ] * 1e6;
%z_valsGPU = [ min(z_distGPU):2.5e-7:max(z_distGPU) ] * 1e6;
[min(y_distGPU) max(y_distGPU) min(z_distGPU) max(z_distGPU)]
y_valsGPU = [ 0:2.5e-7:20.99e-6] * 1e6;
z_valsGPU = [-20.0e-6:2.5e-7:8.1e-6] * 1e6;
%y_valsGPU = [ min(y_distGPU):1e-6:max(y_distGPU) ] * 1e6;
%z_valsGPU = [ min(z_distGPU):1e-6:max(z_distGPU) ] * 1e6;
[min(y_distGPU) max(y_distGPU) min(z_distGPU) max(z_distGPU)]
% force contours at 0.0, 0.1, ..., 0.9+
line_valuesGPU = 0:0.1:0.9;   
[length(p_trapGPU) length(z_valsGPU) length(y_valsGPU) (length(z_valsGPU)*length(y_valsGPU)) length(p_trapGPU)-(length(z_valsGPU)*length(y_valsGPU))]
[B,IX] = sort(y_distGPU,1);
p_trapGPU = p_trapGPU(IX);
% reshape single column of trapping probabilities into a grid
c_matrixGPU = reshape(p_trapGPU(1:length(z_valsGPU)*length(y_valsGPU)), length(z_valsGPU), length(y_valsGPU));

figure,
fontsize = 20;
% plots the contours prescribed by line_values, flipping the Y-axis
xlabel('Y (in \mum)','FontSize',fontsize);
ylabel('Z (in \mum)','FontSize',fontsize);
title('Trapping probability for static trap','FontSize',fontsize);
axis([min(y_valsGPU) max(y_valsGPU) -9.0 max(z_valsGPU)])
%title('Trapping probability for moving trap (0.65\mum/s in Y)','FontSize',fontsize);%, 'Interpreter', 'LaTex')
hold on;
contourf(y_valsGPU, z_valsGPU, c_matrixGPU, line_valuesGPU,'LineStyle',':','LineWidth',1e-100);
set(gca, 'FontSize',fontsize);
set(gca, 'YDir', 'reverse')
hold off;



