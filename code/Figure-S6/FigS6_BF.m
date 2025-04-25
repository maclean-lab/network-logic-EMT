
clear
clc

load param_PertBA0_25pc_OR.mat
load param_SameSA0SB0_OR.mat

low=tri_type_lowBA0;
high=tri_type_highBA0;

nr=length(tri_type_highBA0);

for ii=1:nr

    if tri_type_lowBA0(ii)==1 
        M_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,4);
        H_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,2);
    elseif tri_type_lowBA0(ii)==3
        M_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,4);
        H_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,2);
    elseif tri_type_lowBA0(ii)==2 
        M_state_lowBA0(ii,1)=SNB_lowBA0(ii,1)-SNB_lowBA0(ii,4);
        H_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,2);
    elseif tri_type_lowBA0(ii)==4 
        M_state_lowBA0(ii,1)=SNB_lowBA0(ii,1)-SNB_lowBA0(ii,4);
        H_state_lowBA0(ii,1)=SNB_lowBA0(ii,3)-SNB_lowBA0(ii,2);
    end


    if tri_type_highBA0(ii)==1 
        M_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,4);
        H_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,2);
    elseif tri_type_highBA0(ii)==3
        M_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,4);
        H_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,2);
    elseif tri_type_highBA0(ii)==2
        M_state_highBA0(ii,1)=SNB_highBA0(ii,1)-SNB_highBA0(ii,4);
        H_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,2);
    elseif tri_type_highBA0(ii)==4
        M_state_highBA0(ii,1)=SNB_highBA0(ii,1)-SNB_highBA0(ii,4);
        H_state_highBA0(ii,1)=SNB_highBA0(ii,3)-SNB_highBA0(ii,2);
    end

end

low_types=[sum(low==1),sum(low==2),sum(low==3),sum(low==4)];
low_pc=low_types.*100./sum(low_types);

high_types=[sum(high==1),sum(high==2),sum(high==3),sum(high==4)];
high_pc=high_types.*100./sum(high_types);


low_norm_h = H_state_lowBA0./H_state_lowBA0;
high_norm_h = H_state_highBA0./H_state_lowBA0;

lowx_h_high_gt = find(high_norm_h>1);
lowx_h_high_lt = find(high_norm_h<1);
lowx_h_high_eq = find(high_norm_h==1);
tot_h=length(lowx_h_high_gt)+length(lowx_h_high_lt)+length(lowx_h_high_eq);


figure(10)
subplot(3,3,1)

% Define the vertices of the triangle
x1 = [0, 300, 300];
y1 = [0, 0, 300];

% Plot the vertices

plot(x1, y1, 'k-'); % Plot the lines connecting the vertices


fill(x1, y1, 'b'); % Fill the triangle with blue color
hold on

% Define the vertices of the triangle
x2 = [0, 0, 300];
y2 = [0, 300, 300];


plot(x2, y2, 'k-'); % Plot the lines connecting the vertices


fill(x2, y2, 'g'); % Fill the triangle with blue color
hold on

plot(H_state_lowBA0,H_state_highBA0)
hold on
plot([1,300],[1,300])

hstate_cont=([length(lowx_h_high_gt),length(lowx_h_high_lt),length(lowx_h_high_eq)].*100./tot_h)




SNB_OG = SNB(nrr,:);


diff_OG_high = SNB_OG-SNB_highBA0;

pc_moved_high = diff_OG_high.*100./SNB_OG;

avg_pc_moved_high = mean(pc_moved_high);


std_err_high = std(pc_moved_high) ./ sqrt(size(pc_moved_high, 1));

custom_colors = ['r', 'g', 'b','yellow'];


subplot(3,3,3)
bar(avg_pc_moved_high, 'FaceColor', 'flat');
for i = 1:numel(avg_pc_moved_high)
    c = custom_colors(i);
    barh(i, avg_pc_moved_high(i), 'FaceColor', c);
    hold on;
    %errorbar(i, avg_pc_moved_high(i), std_err_high(i), 'k', 'LineWidth', 1);
    errorbar(avg_pc_moved_high(i), i, std_err_high(i), 'horizontal', 'Color', 'k');
end
hold off;

% Reverse the x-axis
set(gca, 'XDir', 'reverse');








