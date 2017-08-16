clear; clf; clc; format compact
data = dlmread('../dat/traj.dat');
vertex_values = strsplit(fileread('../dat/legend.dat'), '\n');

clf; hold on;

legentries = {};
for i = 2:(size(data,2)-1)
    plot(data(:,1), data(:,i), 'LineWidth', 2)
    legentries{end+1} = [vertex_values{i+1} ': $p_\infty =$ ' sprintf('%0.3f',data(end,i))];
end

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
print(fig, '../img/p_of_t.pdf', '-dpdf')