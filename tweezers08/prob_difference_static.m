%function [] = prob_difference_static(fname1, fname2, title_text, save_fname)

  % files containing [y-distance-to-laser z-distance-to-laser prob-of-trap]
  % These files will be compared to each other using some method (subtraction, whatever) to get an estimate of errors
  % e.g., filtered_cpu vs. static_DOUBLE cpu
%path1 = fname1;
%path2 = fname2;

format long;
clear all;
close all;
clc;
title_text = 'Plot Comparison';
save_fname = 'versus_plot.png';

%path1 = 'outputCPUFull.txt';
%path2 = 'outputGPUFull.txt';
%path1 = 'final1024CPU.txt';
%path2 = 'final1024GPU.txt';
%path1 = 'Results/outputGPU100.txt';
%path2 = 'Results/final1024GPU.txt';

path1 = 'outputCPUFull.txt';
path2 = 'outputCPUDouble.txt';


% Third argument is the title of the plot
if nargin < 3
title_text = 'Plot Comparison'
end

% Fourth argument is the filename of the image to output
if nargin < 4
save_fname = 'versus_plot.png'
end

% loads in lines of three floats, ignoring lines beginning with #
[y_dist1 z_dist1 p_trap1] = textread(path1, '%f %f %f', 'commentstyle', 'shell');
[y_dist2 z_dist2 p_trap2] = textread(path2, '%f %f %f', 'commentstyle', 'shell');

%[min(y_dist1) max(y_dist1) min(z_dist1) max(z_dist1)]
%[min(y_dist2) max(y_dist2) min(z_dist2) max(z_dist2)]

%y_vals = [ 0:2.5e-7:20.99e-6] * 1e6;
%z_vals = [-20.0e-6:2.5e-7:8.1e-6] * 1e6;

% sets the grid size for the contour plot
% if the two grids are different sizes, take the intersection
% and only plot that
max_ymin = max( min(y_dist1), min(y_dist2));
max_zmin = max( min(z_dist1), min(z_dist2));
min_ymax = min( max(y_dist1), max(y_dist2));
min_zmax = min( max(z_dist1), max(z_dist2));

[max_ymin max_zmin min_ymax min_zmax];

y_vals = [ 0:2.5e-7:20.99e-6] * 1e6;
z_vals = [-20.0e-6:2.5e-7:8.1e-6] * 1e6;

%y_vals = [ max_ymin:2.5e-7:min_ymax ] * 1e6;
%z_vals = [ max_zmin:2.5e-7:min_zmax ] * 1e6;

% force contours at 0.0, 0.1, ..., 0.9+
line_values = 0:0.1:0.9;   

% % delete elements outside the intersection of the two bounds
%out_of_bounds1 = find( y_dist1>min_ymax | y_dist1<max_ymin | z_dist1>min_zmax | z_dist1<max_zmin);
%out_of_bounds2 = find( y_dist2>min_ymax | y_dist2<max_ymin | z_dist2>min_zmax | z_dist2<max_zmin);
% p_trap1(out_of_bounds1) = []; y_dist1(out_of_bounds1) = []; z_dist1(out_of_bounds1) = [];
% p_trap2(out_of_bounds2) = []; y_dist2(out_of_bounds2) = []; z_dist2(out_of_bounds2) = [];


[B1,IX1] = sort(y_dist1,1);
p_trap1 = p_trap1(IX1);

[B2,IX2] = sort(y_dist2,1);
p_trap2 = p_trap2(IX2);

% now form the actual trapping values we want to display
% This is some difference between the two sets of trapping values passed in.
p_trap = abs(p_trap2 - p_trap1);

% reshape single column of trapping probabilities into a grid
c_matrix = reshape(p_trap, length(z_vals), length(y_vals));

fontsize = 20;
% plots the contours prescribed by line_values, flipping the Y-axis
xlabel('Y (in \mum)','FontSize',fontsize);
ylabel('Z (in \mum)','FontSize',fontsize);
title(title_text,'FontSize',fontsize);
axis([min(y_vals) max(y_vals) -9.0 max(z_vals)])
%title('Trapping probability for moving trap (0.65\mum/s in Y)','FontSize',fontsize);%, 'Interpreter', 'LaTex')
hold on;
f = contourf(y_vals, z_vals, c_matrix, line_values,'LineStyle',':','LineWidth',1e-100);
set(gca, 'FontSize',fontsize);
set(gca, 'YDir', 'reverse');

%set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6])
if nargin > 3
print('-dpng', save_fname);
end

hold off;

error=abs(p_trap2-p_trap1);

figure, hist(error,200)

%print errors
rmsdError=sqrt((sum(sum(((error).^2))))/(size(error,1)*size(error,2)))
frobeniusError=sqrt(sum(sum(error)))
maxError=max(max(error))
absoluteErrorL1=sum(sum(error))
meanErrorL1=absoluteErrorL1/(size(error,1)*size(error,2))
relativeErrorL1=absoluteErrorL1/sum(sum(p_trap1))

absoluteErrorL2=sqrt(sum(sum(error.^2)))
meanErrorL2=absoluteErrorL2/(size(error,1)*size(error,2))
relativeErrorL2=absoluteErrorL2/sqrt(sum(sum(p_trap1.^2)))



standardDeviation=std(error)


%L1RelativeError=norm(p_trap2-p_trap1,1)/norm(p_trap1,1)
%L2RelativeError=norm(p_trap2-p_trap1,2)/norm(p_trap1,2)

zeros=sum((error)<.000001)
nonZero=(size(error,1)*size(error,2))-zeros
nonZeroPercentage=(nonZero/(size(error,1)*size(error,2)))*100

indices=find(error>.000001);
sum(error(indices))
mean(error(indices))


%A=[1 2 3; 4 5 6; 7 8 10]
%AHat=[.44 2.36 3.04; 3.09 5.87 6.66; 7.36 7.77 9.07]
%norm(AHat-A)/norm(A)
