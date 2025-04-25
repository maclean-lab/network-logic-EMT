
clear
clc

% This code is for the SDE simulations for the AND model
% Run this code for different values of S and store the output as 
% stochEM_sigma1_Unpert_AND.mat, stochEM_sigma1_25pcInd_AND.mat 
% stochEM_sigma1_25pcDir_AND.mat

function [dA, dB] = EMT_AND(A, B, s)
    
    
ksa=87.2400;
ksb=87.2400;
kab=1;
kba=1;
kbb=1;
ab0=400;
ba0=387.1077;
bb0=100;
nsa=1;
nsb=1;
nab=3;
nba=5;
nbb=2;
 ga=0.0773;
 gb=0.0773;
sa0=120;
sb0=120;
ka0=0.9486;
kb0=0.9486;



%sb0 = 120-0.25*120;     % Direct perturbation
    
%sa0 = 120-0.25*120;     % Indirect perturbation


    dA = ka0+(ksa./(1+(s./sa0).^nsa).*kba./(1+(B./ba0).^nba)) - ga.*A;

    dB = kb0+ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb) ...
        .*kbb.*(B./bb0).^nbb./(1+(B./bb0).^nbb) ...
        .*kab./(1+(A./ab0).^nab)-gb.*B;

end



% Solving the above function using Euler-Maruyama method
function [A, B] = euler_maruyama(EMT_AND, A0, B0, s, sigma_A, sigma_B, dt, T)

    N = floor(T / dt);  % Number of time steps
    A = zeros(N, 1);    % Array to store A values
    B = zeros(N, 1);    % Array to store B values
    A(1) = A0;          % Initial condition for A
    B(1) = B0;          % Initial condition for B

    for n = 1:(N - 1)
        dW_A = sqrt(dt) * randn();  % Wiener increment for A (random normal variable)
        dW_B = sqrt(dt) * randn();  % Wiener increment for B (random normal variable
        
        % Get the deterministic parts for both A and B 
        [dA, dB] = EMT_AND(A(n), B(n), s);

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


% Method Parameters
A0 = 500.0;    % Initial value of A
B0 = 0.0;           % Initial value of B          
sigma_A = 1;        % Noise level for A
sigma_B = 1;        % Noise level for B
dt = 0.01;          % Time step
T = 1000.0;         % Total time


sig = [90,100,110,120,130,140,150,160];    % Signal (S) value 

% Loop for the signal value
for ij = 1:length(sig)
    s = sig(ij);

    ij


    % Loop for trajectories at a given signal value

    nr = 1000 ;     % number of trajectories (cells)

    for jj = 1:nr

       
        [A_values, B_values] = euler_maruyama(@EMT_AND, A0, B0, s, sigma_A, sigma_B, dt, T);
        
        A_ss(jj,ij) = A_values(end);
        B_ss(jj,ij) = B_values(end);
        ss_sig(jj,ij) = s;

    end

end

B_sig = B_ss(:);
s_sig = ss_sig(:);
A_sig = A_ss(:);


load Bifur_Unpert_AND.mat
plot(sortedBSSS(:,1),sortedBSSS(:,2)*0.25,'.')
hold on
plot(sortedBUSS(:,1),sortedBUSS(:,2)*0.25,'.')

plot(s_sig,B_sig,'o')

