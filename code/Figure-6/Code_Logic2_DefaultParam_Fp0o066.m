
% This code is for Logic 2 for the parameters from Logic 1 with Fp = 0.066

% Main script
clear all;
close all;

% for WT; vsg1 = 1.202, vsg2 = 1;
% For MT; vsg1 = 1.022, vsg2 = 0.85

% Parameters
Kig = 2; vsg1 = 1.202; vsg2 = 1; Kag1 = 0.28; Kag2 = 0.55; kdg = 1;
Kin2 = 2; vsn1 = 0.856; Kin1 = 0.28; vsn2 = 1; Kan = 0.55; kdn = 1;
vsfr1 = 2.8; Kifr = 0.5; vsfr2 = 2.8; Kafr = 0.5; kdfr = 1;
va = 20; Kd = 2; Ka = 0.7; vi = 3.3; Ki = 0.7;
vex = 0; vsf = 0.6; Kaf = 5; kdf = 0.09;
q = 4; r = 3; s = 4; w = 4; u = 3; v = 4; z = 4;

% For mutant;
%vsg1 = 1.022; vsg2 = 0.85;



% Grid size and variables
n_rows = 5;
n_cols = 5;
n_cells = n_rows * n_cols;
n_vars = 5;
n_loops = 20;

for ij=1:n_loops
    ij

% Gamma values (randomly chosen from uniform distribution [-0.1, 0.1])
rng('shuffle'); % Randomize seed for different runs
gamma = -0.1 + (0.1 - (-0.1)) * rand(n_rows, n_cols); % Uniform distribution [-0.1, 0.1]

% Initial conditions for all 25 cells (G1=0, N1=0 for all)
initial_conditions = zeros(n_cells * n_vars, 1);
for i = 1:n_cells
    idx = (i-1)*n_vars + 1;
    initial_conditions(idx) = 0;      % G1
    initial_conditions(idx+1) = 0;    % N1
    initial_conditions(idx+2) = 2.8;  % FR1
    initial_conditions(idx+3) = 0.25; % ERK1
    initial_conditions(idx+4) = 0.066;% Fs1
end

% Time span
tspan = [0 60];

% Solve ODEs
[t, y] = ode45(@(t,y) odesystem(t,y,gamma,Kig,vsg1,vsg2,Kag1,Kag2,kdg,...
    Kin2,vsn1,Kin1,vsn2,Kan,kdn,vsfr1,Kifr,vsfr2,Kafr,kdfr,va,Kd,Ka,vi,Ki,...
    vex,vsf,Kaf,kdf,q,r,s,w,u,v,z,n_rows,n_cols), tspan, initial_conditions);

% Reshape solution
y_reshaped = reshape(y, [], n_vars, n_cells);

% Create figure with 5x5 subplots
figure('Position', [100, 100, 1000, 800]);
for i = 1:n_cells
    row = ceil(i/n_cols);
    col = mod(i-1, n_cols) + 1;
    subplot(n_rows, n_cols, (row-1)*n_cols + col);
    
    plot(t, y_reshaped(:,1,i), 'b-', 'LineWidth', 1.5); % G1 in blue
    hold on;
    plot(t, y_reshaped(:,2,i), 'r-', 'LineWidth', 1.5); % N1 in red
    hold off;
    
    %title(sprintf('Cell (%d,%d)\n\\gamma=%.3f', row, col, gamma(row,col)));
    % if row == n_rows, xlabel('Time'); end
    % if col == 1, ylabel('Concentration'); end
    % if row == 1 && col == n_cols, legend('G1', 'N1'); end
    % grid on;
    
    % % Verify initial conditions and show gamma
    % fprintf('Cell %d (%d,%d): Initial G1 = %.10f, Initial N1 = %.10f, Gamma = %.6f\n', ...
    %     i, row, col, y_reshaped(1,1,i), y_reshaped(1,2,i), gamma(row,col));
end

% Adjust layout
sgtitle('Logic 2, Default Tristable, WT, Fp = 0.066');
tight_layout();

end

% ODE function
function dydt = odesystem(t,y,gamma,Kig,vsg1,vsg2,Kag1,Kag2,kdg,Kin2,...
    vsn1,Kin1,vsn2,Kan,kdn,vsfr1,Kifr,vsfr2,Kafr,kdfr,va,Kd,Ka,vi,Ki,...
    vex,vsf,Kaf,kdf,q,r,s,w,u,v,z,n_rows,n_cols)

n_cells = n_rows * n_cols;
n_vars = 5;
y = reshape(y, n_vars, n_cells);

% Extract variables
G1 = y(1,:); N1 = y(2,:); FR1 = y(3,:); ERK1 = y(4,:); Fs1 = y(5,:);

% Reshape Fs1 into 5x5 grid for neighbor calculations
Fs1_grid = reshape(Fs1, n_rows, n_cols);

% Calculate F for each cell based on nearest neighbors
F = zeros(n_rows, n_cols);
for i = 1:n_rows
    for j = 1:n_cols
        neighbors = [];
        % Check top neighbor
        if i > 1, neighbors = [neighbors, Fs1_grid(i-1,j)]; end
        % Check bottom neighbor
        if i < n_rows, neighbors = [neighbors, Fs1_grid(i+1,j)]; end
        % Check left neighbor
        if j > 1, neighbors = [neighbors, Fs1_grid(i,j-1)]; end
        % Check right neighbor
        if j < n_cols, neighbors = [neighbors, Fs1_grid(i,j+1)]; end
        % Average of neighbors (if no neighbors, use own Fs1)
        if ~isempty(neighbors)
            F(i,j) = mean(neighbors);
        else
            F(i,j) = Fs1_grid(i,j);
        end
    end
end

% Calculate Fpi for each cell
Fp = (1 + gamma) .* F;
Fp = Fp(:)';  % Convert to row vector matching cell indexing

% Initialize derivatives
dydt = zeros(n_vars, n_cells);

% ODEs for each cell
for i = 1:n_cells
    % dG1/dt
    dydt(1,i) = vsg1*ERK1(i)^r/(Kag1^r+ERK1(i)^r)+((vsg2*G1(i)^s/(Kag2^s+G1(i)^s))*(Kig^q/(Kig^q+N1(i)^q)))-kdg*G1(i);
    
    % dN1/dt
    dydt(2,i) = vsn1*Kin1^u/(Kin1^u+ERK1(i)^u)+((vsn2*N1(i)^v/(Kan^v+N1(i)^v))*(Kin2^w/(Kin2^w+G1(i)^w)))-kdn*N1(i);
    
    % dFR1/dt
    dydt(3,i) = vsfr1*Kifr/(Kifr+N1(i))+(vsfr2*G1(i))/(Kafr+G1(i))-kdfr*FR1(i);
    
    % dERK1/dt
    dydt(4,i) = va*FR1(i)*Fp(i)/(Kd+Fp(i))*((1-ERK1(i))/(Ka+1-ERK1(i)))-vi*ERK1(i)/(Ki+ERK1(i));
    
    % dFs1/dt
    dydt(5,i) = vex+(vsf*N1(i)^z)/(Kaf^z+N1(i)^z)-kdf*Fs1(i);
end

% Reshape output
dydt = dydt(:);
end


% Helper function for tight layout
function tight_layout()
    if ~exist('tightfig', 'file')
        set(gcf,'Position',[100 100 1000 800]);
        axes = findobj(gcf,'Type','axes');
        for ax = axes'
            outerpos = get(ax,'OuterPosition');
            set(ax,'Position',[outerpos(1)+0.02 outerpos(2)+0.02 ...
                outerpos(3)-0.04 outerpos(4)-0.04]);
        end
    else
        tightfig;
    end
end

