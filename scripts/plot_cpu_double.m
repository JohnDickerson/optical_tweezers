% Takes the cpu_double_results file---as one big file, not the 28-odd
% separate ones generated by a cluster run---and plots a contour/heatmap of
% the probabilities across y/z initial particle positions

% Note that the files are [<y> <x> <probability>], with #-delim comments
[y z prob] = textread('cpu_double_results.txt', '%f %f %f', 'commentstyle', 'shell');

% Sort the data by y-values first (col 1), then z-values (col 2) flipped
ordered = sortrows([y z prob], [1 -2]);

% With a step of 0.25 and [0,21) y-values, [-20,8] z-values, we need to
% have an 84 x 113 matrix for the contour plot
y_mat = reshape(ordered(:,1), 113, 84);
z_mat = reshape(ordered(:,2), 113, 84);
prob_mat = reshape(ordered(:,3), 113, 84);

% Plot the matrix as a contour map with N levels; red is 1, blue is 0 prob.
num_levels = 9;

%hold on
contourf(y_mat, -z_mat, prob_mat, num_levels, ':')
axis equal;
axis([0*1e-6 22*1e-6 -8*1e-6 10*1e-6 ]);
font_size = 14;
xlabel('Y', 'FontSize', font_size);
ylabel('Z', 'FontSize', font_size);
title('Trapping Probabilities for Static Trap, CPU Double Precision', 'FontSize', font_size);
%hold off