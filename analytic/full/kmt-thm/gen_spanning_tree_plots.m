clear; clc;

L = 6;

A = textread(['adjmx_' num2str(L) '.dat']);
v = textread(['vtx_' num2str(L) '.dat']);
nv = size(v,1);
A = sparse(A);


[idx, src, dst] = lib.generateSpanningTrees(A);

treefile = fopen(['trees_' num2str(L) '.dat'], 'w');

disp(size(idx,2))
for i = 1:size(idx, 2)

	f = figure('visible', 'off');
	hold on;
	gplot(A, v, '-ob')
	axis square; grid on; box on;
	set(gca, 'FontSize', 16)
	
    for j = 1:size(idx, 1)
        
		fprintf(treefile, '%i %i', src(idx(j, i)), dst(idx(j, i)));
		
		coords = [v(src(idx(j, i)),:); v(dst(idx(j, i)),:)];
		x_coords = coords(:,1);
		y_coords = coords(:,2);
		
		plot(x_coords, y_coords, 'r', 'LineWidth', 4)
		
		if j ~= size(idx, 2)
			fprintf(treefile, ' ');
		end
    end

    fprintf(treefile, '\n');
	f.PaperPositionMode = 'auto';
	fig_pos = f.PaperPosition;
	f.PaperSize = [fig_pos(3) fig_pos(4)];
	print(f, ['tree_' num2str(L) '_' num2str(i) '.pdf'], '-dpdf')
	close(f)

end

fprintf('\n')
