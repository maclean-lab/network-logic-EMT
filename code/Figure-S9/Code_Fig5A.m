
clear
clc

% This MATLAB code is to plot Figure 5A

load stochEM_sigma1_Unpert_OR.mat

res_aa = A_sig;
res_bb = B_sig;
res_rr = s_sig;

max_aa = max(res_aa);
min_aa = min(res_aa);
max_bb = max(res_bb);
min_bb = min(res_bb);

res_aa1 = (res_aa-min_aa)./(max_bb-min_aa);
res_bb1 = (res_bb-min_bb)./(max_bb-min_bb);

subplot(2,3,1)
plot(res_rr,res_bb1,'.')
hold on
