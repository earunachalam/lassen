%% calculate the distribution of segment lengths for different alpha
% where alpha is the driven rate constant
% for systems with L = 20 (total number of monomers)

clear; clc; format compact
clf; hold on;

p = []; f = {}; g = {};

alphavals = {'0e0', '1e-1', '1e0', '2e0', '5e0'};
for alpha = alphavals
    
    data = dlmread(['../dat/traj_a_' alpha{1} '_L_20.dat']);
    ssval = data(end,2:end); % discard t column

    vertex_values = strsplit(fileread('../dat/legend_L_20.dat'), '\n');
    vertex_values(end) = [];
    vertex_values_numeric = vertex_values;
    vertex_length_avg = nan(1, numel(vertex_values));
    pl = zeros(1, str2double(vertex_values{end}));

    for i = 1:numel(vertex_values)

        strsplit(vertex_values{i}, ',');
        vertex_values_numeric{i} = strsplit(vertex_values{i}, ',');

        for j = 1:numel(vertex_values_numeric{i})
            vertex_values_numeric{i}{j} = str2double(vertex_values_numeric{i}{j});
        end

        curr_vertex_value = cell2mat(vertex_values_numeric{i});
        vertex_values_numeric{i} = curr_vertex_value;
        vertex_length_avg(i) = mean(vertex_values_numeric{i});

        num_total_segs = numel(curr_vertex_value);
        for j = 1:numel(pl)
            num_matching_segs = sum(curr_vertex_value == j);
            p_matching_segs = num_matching_segs/num_total_segs;
            pl(j) = pl(j) + p_matching_segs*ssval(i);
        end
    end

    l = 1:length(pl);
    p(end+1) = plot(l, pl, 'LineWidth', 2, 'DisplayName', ['$\alpha =$ ' num2str(str2double(alpha{1})) ', $\bar{l} = $' sprintf('%0.3f', sum(vertex_length_avg .* ssval))]);
    [f{end+1}, gof] = fit(l(2:end-1)',pl(2:end-1)','exp1');
%     [g{end+1}, gof2] = fit(l(2:end-1)',pl(2:end-1)','power1')
    [g{end+1}, gof2] = fit(log(l(2:end-1)'),log(pl(2:end-1)'),'poly1')
%     [g{end+1}, gof2] = fit(log(l(2:end-1)'),log(pl(2:end-1)'), @(a,b,x) a*b.^(-x), 'StartPoint',[0 0])
    domain_fit = linspace(2, length(pl)-1, 1000);
    cv = coeffvalues(g{end}); plot(domain_fit, exp(cv(2) + domain_fit.*cv(1)), '--', 'LineWidth', 2)
%     plot(g{end})
    return
end

set(gca, 'Fontsize', 16)
box on
xlabel('segment length $l$', 'Interpreter', 'latex')
ylabel('$p(l)$', 'Interpreter', 'latex')
leg = legend(p);
set(leg, 'Interpreter', 'latex', 'Location', 'East')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/p_of_l.pdf', '-dpdf')

alphavals_copy = alphavals;

%% get the values of the exponential decay fit coefficients

alphavals = alphavals_copy;
a = nan(size(f));
b = a;
for i = 1:numel(f)
    c = coeffvalues(f{i});
    a(i) = c(1);
    b(i) = c(2);
    alphavals{i} = str2double(alphavals_copy{i});
end

alphavals = cell2mat(alphavals);

%% plot pre-exponential factor a as a function of driven rate constant alpha

clf
hold on
scatter(alphavals, a, 54, 'filled')
set(gca, 'Fontsize', 16)
box on
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$a$ (in $a \exp(b \alpha) $)', 'Interpreter', 'latex')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/preexp_vs_alpha.pdf', '-dpdf')

%% plot exponential argument b as a function of driven rate constant alpha

clf
hold on
scatter(alphavals, b, 54, 'filled')
set(gca, 'Fontsize', 16)
box on
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$b$ (in $a \exp(b \alpha) $)', 'Interpreter', 'latex')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/exparg_vs_alpha.pdf', '-dpdf')

%% get the values of the power law fit coefficients

alphavals = alphavals_copy;
a = nan(size(g));
b = a;
for i = 1:numel(g)
    c = coeffvalues(g{i});
    a(i) = c(1);
    b(i) = c(2);
    alphavals{i} = str2double(alphavals_copy{i});
end

alphavals = cell2mat(alphavals);

%% plot prefactor a as a function of driven rate constant alpha

clf
hold on
scatter(alphavals, a, 54, 'filled')
set(gca, 'Fontsize', 16)
box on
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$a$ (in $a \alpha^b $)', 'Interpreter', 'latex')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/prefactor_vs_alpha.pdf', '-dpdf')

%% plot exponent b as a function of driven rate constant alpha

clf
hold on
scatter(alphavals, b, 54, 'filled')
set(gca, 'Fontsize', 16)
box on
xlabel('$\alpha$', 'Interpreter', 'latex')
ylabel('$b$ (in $a \alpha^b $)', 'Interpreter', 'latex')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../img/exp_vs_alpha.pdf', '-dpdf')