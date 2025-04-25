
clear
clc

% This code is for the SDE simulations for the OR model
% Run this code for different values of S and store the output as 
% stochEM_sigma1_Unpert_OR.mat, stochEM_sigma1_25pcInd_OR.mat 
% stochEM_sigma1_25pcDir_OR.mat


function [dA, dB] = coupled_system(A, B, s)
    
    
    ksa=8;
    ksb=8;
    kab=10;
    kba=10;
    kbb=7;
    sa0=230;
    sb0=230;
    ab0=80;
    ba0=80;
    bb0=36;
    nsa=4;
    nsb=4;
    nab=6;
    nba=6;
    nbb=3;
    ga=0.13;
    gb=0.13;

    %sb0 = 230-0.25*230;    % Direct perturbation
    %sa0 = 230-0.25*230;    % Indirect perturbation


    dA = (ksa./(1+(s./sa0).^nsa)+kba./(1+(B./ba0).^nba)) - ga.*A;

    dB = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb) ...
        +kbb.*(B./bb0).^nbb./(1+(B./bb0).^nbb) ...
        +kab./(1+(A./ab0).^nab)-gb.*B;

end

% Euler-Maruyama method 
function [A, B] = euler_maruyama(coupled_system, A0, B0, s, sigma_A, sigma_B, dt, T)
    N = floor(T / dt);  % Number of time steps
    A = zeros(N, 1);    % Array to store A values
    B = zeros(N, 1);    % Array to store B values
    A(1) = A0;          % Initial condition for A
    B(1) = B0;          % Initial condition for B

    for n = 1:(N - 1)
        dW_A = sqrt(dt) * randn();  % Wiener increment for A (random normal variable)
        dW_B = sqrt(dt) * randn();  % Wiener increment for B (random normal variable)
        
        
        % Get the deterministic parts for both A and B 
        [dA, dB] = coupled_system(A(n), B(n), s);

        if A(n)<0
            A(n) = 0;
        end

        if B(n)<0
            B(n) = 0;
        end
        
        % Update A and B using Euler-Maruyama method
        A(n + 1) = A(n) + dA * dt + sigma_A * dW_A;
        B(n + 1) = B(n) + dB * dt + sigma_B * dW_B;
    end
end


% Parameters
A0 = 500.0;         % Initial value of A
B0 = 0.0;           % Initial value of B
%s = 10;           
sigma_A = 1;        % Noise level for A
sigma_B = 1;        % Noise level for B
dt = 0.01;          % Time step
T = 1000.0;         % Total time


% Signal (S) value for which SDE is performed
sig = [90,100,110,120,130,140,150,160];

% Loop for the signal value
for ij = 1:length(sig)
    s = sig(ij);

    % Loop for trajectories at a given signal value

    nr = 1000 ;     % number of trajectories (cells)

    for jj = 1:nr

        
        [A_values, B_values] = euler_maruyama(@coupled_system, A0, B0, s, sigma_A, sigma_B, dt, T);
        
        A_ss(jj,ij) = A_values(end);
        B_ss(jj,ij) = B_values(end);
        ss_sig(jj,ij) = s;

    end

end

B_sig = B_ss(:);
s_sig = ss_sig(:);
A_sig = A_ss(:);


load BifurPlot_Unpert_OR.mat
plot(bdia(:,1),bdia(:,2))
hold on
plot(s_sig,B_sig,'o')



