
% This code plots the shift in the SN points in the perturbed models
% relative to the unperturbed models. The SN points for the perturbed
% models are stored in 'param_ThrePert_SA0SB0_25pc_AND.mat' and the SN
% points for the unperturbed models are stored in
% 'param_RandSA0SB0_AND.mat'.

clear
clc

load param_ThrePert_SA0SB0_25pc_AND.mat
load param_RandSA0SB0_AND.mat

SNB_OG = SNB(nrr,:);    % SN points for unperturbed models

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of the SN points for the direct perturbations

diff_OG_dir = SNB_OG-SNB_SB0ltSA0;

pc_moved_dir = diff_OG_dir.*100./SNB_OG;

avg_pc_moved_dir = mean(pc_moved_dir);


std_err_dir = std(pc_moved_dir) ./ sqrt(size(pc_moved_dir, 1));

custom_colors = ['r', 'g', 'b','yellow'];

figure(2)
subplot(3,2,5)
bar(avg_pc_moved_dir, 'FaceColor', 'flat');

for i = 1:numel(avg_pc_moved_dir)
    c = custom_colors(i);
    barh(i, avg_pc_moved_dir(i), 'FaceColor', c);
    hold on;
    errorbar(avg_pc_moved_dir(i), i, std_err_dir(i), 'horizontal', 'Color', 'k');
end

hold off;


% Reverse the x-axis
set(gca, 'XDir', 'reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of the SN points for the indirect perturbations

diff_OG_ind = SNB_OG-SNB_SA0ltSB0;

pc_moved_ind = diff_OG_ind.*100./SNB_OG;

avg_pc_moved_ind = mean(pc_moved_ind);


std_err_ind = std(pc_moved_ind) ./ sqrt(size(pc_moved_ind, 1));

subplot(3,2,6)
bar(avg_pc_moved_ind, 'FaceColor', 'flat');

for i = 1:numel(avg_pc_moved_ind)
    c = custom_colors(i);
    barh(i, avg_pc_moved_ind(i), 'FaceColor', c);
    hold on;
    errorbar(avg_pc_moved_ind(i), i, std_err_ind(i), 'horizontal', 'Color', 'k');
end

hold off;

% Reverse the x-axis
set(gca, 'XDir', 'reverse');