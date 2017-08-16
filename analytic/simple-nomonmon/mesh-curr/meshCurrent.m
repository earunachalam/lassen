clear; clc; format compact;

syms a b

% initialize arrays with dummy symbolic expression to make sure arrays hold
% syms rather than doubles
eqs = [0 == a+b];   % array to hold all loop equations
all_I = [a];        % array to hold all loop currents

Ncols = 5;

horiz_range = 1:Ncols;

for coord_horiz = 1:Ncols
    coord_vert_max = ceil(coord_horiz/2);
    for coord_vert = 1:coord_vert_max
        vtx_locn = [num2str(coord_horiz) '_' num2str(coord_vert)];
        eval(['syms C_' vtx_locn]); % create central loop variable (big I coeff)
        eval(['syms H_' vtx_locn]); % create neighbor loop variable (small I coeff)
        eval(['syms I_' vtx_locn]); % create current variable (I)
        eval(['all_I(end+1) = I_' vtx_locn ';']);
    end
end

% construct loop equations
for coord_horiz = 1:Ncols
    coord_vert_max = ceil(coord_horiz/2);
    for coord_vert = 1:coord_vert_max
        
        left    = coord_horiz - 1;
        right   = coord_horiz + 1;
        up      = coord_vert - 1;
        down    = coord_vert + 1;
        
        vtx_locn = [num2str(coord_horiz) '_' num2str(coord_vert)];
        eval(['eqs(end+1) = 0 == C_' vtx_locn '*I_' vtx_locn ';']);
        
        % add terms from western neighbor
        west_neighbor   = ['H_' num2str(left) '_' num2str(coord_vert)];
        west_current    = ['I_' num2str(left) '_' num2str(coord_vert)];
        if exist(west_neighbor, 'var')
            eval(['eqs(end) = 0 == rhs(eqs(end)) + ' west_neighbor '*' west_current ';']);
        end
        
        % add terms from eastern neighbor
        east_neighbor   = ['H_' num2str(right) '_' num2str(coord_vert)];
        east_current    = ['I_' num2str(right) '_' num2str(coord_vert)];
        if exist(east_neighbor, 'var')
            eval(['eqs(end) = 0 == rhs(eqs(end)) + ' east_neighbor '*' east_current ';']);
        end
        
        % add terms from southern neighbor
        south_neighbor   = ['H_' num2str(coord_horiz) '_' num2str(down)];
        south_current    = ['I_' num2str(coord_horiz) '_' num2str(down)];
        if exist(south_neighbor, 'var')
            eval(['eqs(end) = 0 == rhs(eqs(end)) + ' south_neighbor '*' south_current ';']);
        end
        
        % add terms from northern neighbor
        north_neighbor   = ['H_' num2str(coord_horiz) '_' num2str(up)];
        north_current    = ['I_' num2str(coord_horiz) '_' num2str(up)];
        if exist(north_neighbor, 'var')
            eval(['eqs(end) = 0 == rhs(eqs(end)) + ' north_neighbor '*' north_current ';']);
        end
    end
end

% remove dummy symbolic variables
eqs(1) = [];
all_I(1) = [];

soln = solve(eqs, all_I)

[A,b] = equationsToMatrix(eqs, all_I);
null(A)