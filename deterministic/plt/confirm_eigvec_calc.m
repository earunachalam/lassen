% confirm calculation for case with no multiplicity, performed in
% Mathematica using two different methods: (1) solving the rate equations
% in MML's formalism explicitly, and (2) constructing the rate matrix and
% diagonalizing, getting eigenvalues and eigenvectors.etc.

% plugged in alpha = 1 for both cases, which give the same result.

clear; clf; clc; format compact
data = dlmread('../dat/traj_nomultiplicity.dat');
vertex_values = strsplit(fileread('../dat/legend_L_4.dat'), '\n');

data_less_one = data(:,3:6);
each_t_sum_p = sum(data_less_one,2);
rep_each_t_sum_p = repmat(each_t_sum_p,1,4);
renormed_p = data_less_one./rep_each_t_sum_p;

clf; hold on;
% axis square
legentries = {};
for i = 1:4
    plot(data(:,1), renormed_p(:,i), 'LineWidth', 2)
    legentries{end+1} = [vertex_values{i+1} ': $p_\infty =$ ' sprintf('%0.3f',renormed_p(end,i))];
end
% renormed_p(end,:)

set(gca, 'Fontsize', 16)
box on
xlabel('$t$/s', 'Interpreter', 'latex')
ylabel('$p(t)$', 'Interpreter', 'latex')
legend(legentries, 'Interpreter', 'latex', 'Location', 'East')
xlim([0 4])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/confirm_nomultiplicity_result.pdf', '-dpdf')