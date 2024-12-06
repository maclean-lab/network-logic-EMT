
clear
clc

t0 = 0;        % Initial time
tfinal = 10000;

% The following are the initial conditions
% For each condition, run these initial conditions one by one

%c0 = [10000; 0; 50];     % mZ lower

%c0 = [0;10000;0];      % mZ upper

%c0 = [600; 600; 600];  % mZ intermediate




sig=[15000:500:250000];



for ik=1:length(sig)
    ik

    s=sig(ik);

% Define ODE function
ode_func = @(t, c) MJ_Model(t,c,s);
% Integrate the ODEs
[t, c] = ode45(ode_func, [t0, tfinal], c0);

% Extract the steady state values of c
ss_mi200(ik,1) = c(end, 1);
ss_mz(ik,1) = c(end,2);
ss_z(ik,1) = c(end,3);

clear t c

end



%subplot(2,2,1)%
figure(3)
subplot(2,2,1)
plot(sig,ss_mz,'.');
hold on

subplot(2,2,3)
plot(sig,ss_z,'.');
hold on
% 


function dcdt = MJ_Model(t,c,s)

dcdt=zeros(3,1);

%s = 250000;

mu0=10000;

gmi = 2100;
gmz = 11;
gz = 100;
kmi = 0.05;
kmz = 0.5;
kz = 0.1;
z0_mi = 220000;
z0_mz = 25000;
s0_mi = 180000;
s0_mz = 180000;
nz_mi = 3;
ns_mi = 2;
nz_mz = 2;
ns_mz = 2;
lambda_z_mi = 0.1;
lambda_s_mi = 0.1;
lambda_z_mz = 7.5;
lambda_s_mz = 10;

s0_mi = 80000;      % for Indirect > Direct

%s0_mz = 80000;      % FOR Direct > Indirect


sig_mz=[0.04,0.2,1,1,1,1];
Ym1=0.0;
for i1=1:6
Ym1=Ym1+sig_mz(i1).*nchoosek(6,i1).*(c(1)./mu0).^i1;
end

Ymz=(1./(1+c(1)./mu0).^6).*Ym1;

sig_mi=[0.005,0.05,0.5,0.5,0.5,0.5];
Ymu1=0;
for i2=1:6
  Ymu1=Ymu1+sig_mi(i2).*i2.*nchoosek(6,i2).*(c(1)./mu0).^i2;      
end

Ymi=(1./(1+c(1)./mu0).^6).*Ymu1;

sig_li=[0.6,0.3,0.1,0.05,0.05,0.05];
L1=1.0;
for i3=1:6
L1=L1+sig_li(i3).*nchoosek(6,i3).*(c(1)./mu0).^i3;
end

L=(1./(1+c(1)./mu0).^6).*L1;


H_z_mi = (1+lambda_z_mi.*(c(3)./z0_mi).^nz_mi)./(1+(c(3)./z0_mi).^nz_mi);

H_s_mi = (1+lambda_s_mi.*(s./s0_mi).^ns_mi)./((1+(s./s0_mi).^ns_mi));

H_z_mz = (1+lambda_z_mz.*(c(3)./z0_mz).^nz_mz)./(1+(c(3)./z0_mz).^nz_mz);

H_s_mz = (1+lambda_s_mz.*(s./s0_mz).^ns_mz)./(1+(s./s0_mz).^ns_mz);


dcdt(1) = gmi.*H_z_mi.*H_s_mi-c(2).*Ymi-kmi.*c(1);

dcdt(2) = gmz.*H_z_mz.*H_s_mz-c(2).*Ymz-kmz.*c(2);

dcdt(3) = gz.*c(2).*L-kz.*c(3);

end