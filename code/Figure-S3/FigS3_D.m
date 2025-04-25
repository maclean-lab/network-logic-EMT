clear
clc

load param_ThrePert_SA0SB0_25pc_AND.mat
load param_RandSA0SB0_AND.mat

param_unp = param_RandSA0SB0(:,nrr);
%param_thre = param_ThrePert_OR(:,nrr);


unp=tri_types(nrr,:);
ind=tri_type_SA0ltSB0;
dir=tri_type_SB0ltSA0;

[unp_sorted,idx] = sort(unp);
dir_sorted = dir(idx);
ind_sorted = ind(idx);

SN_DirSorted = SNB_SB0ltSA0(idx,:);
SN_IndSorted = SNB_SA0ltSB0(idx,:);

SNB_Unp = SNB(nrr,:);
SNB_Dir = SNB_SB0ltSA0;
SNB_Ind = SNB_SA0ltSB0;

% Type 4 (1U1D): SN1>SN3>SN2>SN4
% Type 2 (1U2D): SN1>SN3>SN4>SN2
% Type 3 (2U1D): SN3>SN1>SN2>SN4
% Type 1 (2U2D): SN3>SN1>SN4>SN2

nr = length(nrr);

for ik = 1:nr

    if SNB_Unp(ik,1)>SNB_Unp(ik,3) && SNB_Unp(ik,3)>SNB_Unp(ik,2) && SNB_Unp(ik,2)>SNB_Unp(ik,4)
        triType_Unp(ik,1) = 4;
    elseif SNB_Unp(ik,1)>SNB_Unp(ik,3) && SNB_Unp(ik,3)>SNB_Unp(ik,4) && SNB_Unp(ik,4)>SNB_Unp(ik,2)
        triType_Unp(ik,1) = 2;
    elseif SNB_Unp(ik,3)>SNB_Unp(ik,1) && SNB_Unp(ik,1)>SNB_Unp(ik,2) && SNB_Unp(ik,2)>SNB_Unp(ik,4)
        triType_Unp(ik,1) = 3;
    elseif SNB_Unp(ik,3)>SNB_Unp(ik,1) && SNB_Unp(ik,1)>SNB_Unp(ik,4) && SNB_Unp(ik,4)>SNB_Unp(ik,2)
        triType_Unp(ik,1) = 1;
    end

    if SNB_Dir(ik,1)>SNB_Dir(ik,3) && SNB_Dir(ik,3)>SNB_Dir(ik,2) && SNB_Dir(ik,2)>SNB_Dir(ik,4)
        triType_Dir(ik,1) = 4;
    elseif SNB_Dir(ik,1)>SNB_Dir(ik,3) && SNB_Dir(ik,3)>SNB_Dir(ik,4) && SNB_Dir(ik,4)>SNB_Dir(ik,2)
        triType_Dir(ik,1) = 2;
    elseif SNB_Dir(ik,3)>SNB_Dir(ik,1) && SNB_Dir(ik,1)>SNB_Dir(ik,2) && SNB_Dir(ik,2)>SNB_Dir(ik,4)
        triType_Dir(ik,1) = 3;
    elseif SNB_Dir(ik,3)>SNB_Dir(ik,1) && SNB_Dir(ik,1)>SNB_Dir(ik,4) && SNB_Dir(ik,4)>SNB_Dir(ik,2)
        triType_Dir(ik,1) = 1;
    end

    if SNB_Ind(ik,1)>SNB_Ind(ik,3) && SNB_Ind(ik,3)>SNB_Ind(ik,2) && SNB_Ind(ik,2)>SNB_Ind(ik,4)
        triType_Ind(ik,1) = 4;
    elseif SNB_Ind(ik,1)>SNB_Ind(ik,3) && SNB_Ind(ik,3)>SNB_Ind(ik,4) && SNB_Ind(ik,4)>SNB_Ind(ik,2)
        triType_Ind(ik,1) = 2;
    elseif SNB_Ind(ik,3)>SNB_Ind(ik,1) && SNB_Ind(ik,1)>SNB_Ind(ik,2) && SNB_Ind(ik,2)>SNB_Ind(ik,4)
        triType_Ind(ik,1) = 3;
    elseif SNB_Ind(ik,3)>SNB_Ind(ik,1) && SNB_Ind(ik,1)>SNB_Ind(ik,4) && SNB_Ind(ik,4)>SNB_Ind(ik,2)
        triType_Ind(ik,1) = 1;
    end


end

diff_ind_SN31 = SNB_Ind(:,3)-SNB_Ind(:,1);
diff_ind_SN42 = SNB_Ind(:,4)-SNB_Ind(:,2);
diff_ind = diff_ind_SN31+diff_ind_SN42;

diff_dir_SN31 = SNB_Dir(:,3)-SNB_Dir(:,1);
diff_dir_SN42 = SNB_Dir(:,4)-SNB_Dir(:,2);
diff_dir = diff_dir_SN31+diff_dir_SN42;

percent_ind = (length(find(diff_ind>diff_dir))/length(diff_ind))*100



percent_ind_EMT = (length(find(diff_ind_SN31>diff_dir_SN31))/length(diff_ind_SN31))*100
percent_ind_MET = (length(find(diff_ind_SN42>diff_dir_SN42))/length(diff_ind_SN42))*100

subplot(4,4,1)
plot(diff_ind_SN31,diff_dir_SN31,'o')
hold on
plot([-250,250],[-250,250])
plot([-250,250],[0,0])
plot([0,0],[-250,250])
xlabel('Ind (SN3-SN1)')
ylabel('Dir (SN3-SN1)')
title('During EMT')


subplot(4,4,2)
plot(diff_ind_SN42,diff_dir_SN42,'o')
hold on
plot([-250,250],[-250,250])
plot([-250,250],[0,0])
plot([0,0],[-250,250])
xlabel('Ind (SN4-SN2)')
ylabel('Dir (SN4-SN2)')
title('During MET')


% Convert data into a matrix
data = [unp_sorted';dir_sorted';ind_sorted'];


% Create a heatmap with default colormap
%heatmap(data, 'ColorbarVisible', 'off');
% figure(22)
% subplot(3,1,1)
% heatmap(data,'colormap',parula);
% xlabel('Models');
% %ylabel('a OR b');
% title('OR models (25% change in S_A_0, S_B_0 thresholds)');
% 
% len_H_Dir = SNB_SB0ltSA0(:,3)-SNB_SB0ltSA0(:,2);
% len_H_Ind = SNB_SA0ltSB0(:,3)-SNB_SA0ltSB0(:,2);
% 
% subplot(3,1,2)
% plot(len_H_Ind,len_H_Dir,'o')
% hold on
% plot([0,200],[0,200])
% xlabel('Indirect')
% ylabel('Direct')
