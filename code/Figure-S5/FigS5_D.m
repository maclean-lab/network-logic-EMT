clear
clc

load param_pertBB0_25pc_OR.mat
load param_SameSA0SB0_OR.mat

param_unp = param_OR(:,nrr);
%param_thre = param_ThrePert_OR(:,nrr);


unp=tri_types(nrr,:);
high=tri_type_highBB0;
low=tri_type_lowBB0;

[unp_sorted,idx] = sort(unp);
low_sorted = low(idx);
high_sorted = high(idx);

SN_lowSorted = SNB_lowBB0(idx,:);
SN_highSorted = SNB_highBB0(idx,:);

SNB_Unp = SNB(nrr,:);
SNB_low = SNB_lowBB0;
SNB_high = SNB_highBB0;

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

    if SNB_low(ik,1)>SNB_low(ik,3) && SNB_low(ik,3)>SNB_low(ik,2) && SNB_low(ik,2)>SNB_low(ik,4)
        triType_low(ik,1) = 4;
    elseif SNB_low(ik,1)>SNB_low(ik,3) && SNB_low(ik,3)>SNB_low(ik,4) && SNB_low(ik,4)>SNB_low(ik,2)
        triType_low(ik,1) = 2;
    elseif SNB_low(ik,3)>SNB_low(ik,1) && SNB_low(ik,1)>SNB_low(ik,2) && SNB_low(ik,2)>SNB_low(ik,4)
        triType_low(ik,1) = 3;
    elseif SNB_low(ik,3)>SNB_low(ik,1) && SNB_low(ik,1)>SNB_low(ik,4) && SNB_low(ik,4)>SNB_low(ik,2)
        triType_low(ik,1) = 1;
    end

    if SNB_high(ik,1)>SNB_high(ik,3) && SNB_high(ik,3)>SNB_high(ik,2) && SNB_high(ik,2)>SNB_high(ik,4)
        triType_high(ik,1) = 4;
    elseif SNB_high(ik,1)>SNB_high(ik,3) && SNB_high(ik,3)>SNB_high(ik,4) && SNB_high(ik,4)>SNB_high(ik,2)
        triType_high(ik,1) = 2;
    elseif SNB_high(ik,3)>SNB_high(ik,1) && SNB_high(ik,1)>SNB_high(ik,2) && SNB_high(ik,2)>SNB_high(ik,4)
        triType_high(ik,1) = 3;
    elseif SNB_high(ik,3)>SNB_high(ik,1) && SNB_high(ik,1)>SNB_high(ik,4) && SNB_high(ik,4)>SNB_high(ik,2)
        triType_high(ik,1) = 1;
    end


end

diff_high_SN31 = SNB_high(:,3)-SNB_high(:,1);
diff_high_SN42 = SNB_high(:,4)-SNB_high(:,2);
diff_high = diff_high_SN31+diff_high_SN42;

diff_low_SN31 = SNB_low(:,3)-SNB_low(:,1);
diff_low_SN42 = SNB_low(:,4)-SNB_low(:,2);
diff_low = diff_low_SN31+diff_low_SN42;

% percent_high = (length(high(diff_high>diff_low))/length(diff_high))*100
% 
% subplot(4,4,7)
% plot(diff_high,diff_low,'o')
% hold on
% plot([-250,250],[-250,250])
% plot([-250,250],[0,0])
% plot([0,0],[-250,250])
% xlabel('high')
% ylabel('low')


percent_high_EMT = (length(high(diff_high_SN31>diff_low_SN31))/length(diff_high_SN31))*100
percent_high_MET = (length(high(diff_high_SN42>diff_low_SN42))/length(diff_high_SN42))*100

subplot(4,4,3)
plot(diff_high_SN31,diff_low_SN31,'o')
hold on
plot([-250,250],[-250,250])
plot([-250,250],[0,0])
plot([0,0],[-250,250])
xlabel('high SA (SN3-SN1)')
ylabel('low SA (SN3-SN1)')
title('During EMT')



subplot(4,4,4)
plot(diff_high_SN42,diff_low_SN42,'o')
hold on
plot([-250,250],[-250,250])
plot([-250,250],[0,0])
plot([0,0],[-250,250])
xlabel('high SA (SN4-SN2)')
ylabel('low SA (SN4-SN2)')
title('During MET')



% Convert data into a matrix
data = [unp_sorted';low_sorted';high_sorted'];


% Create a heatmap with default colormap
%heatmap(data, 'ColorbarVisible', 'off');
% figure(22)
% subplot(3,1,1)
% heatmap(data,'colormap',parula);
% xlabel('Models');
% %ylabel('a OR b');
% title('OR models (25% change in S_A_0, S_B_0 thresholds)');
% 
% len_H_low = SNB_SB0ltSA0(:,3)-SNB_SB0ltSA0(:,2);
% len_H_high = SNB_SA0ltSB0(:,3)-SNB_SA0ltSB0(:,2);
% 
% subplot(3,1,2)
% plot(len_H_high,len_H_low,'o')
% hold on
% plot([0,200],[0,200])
% xlabel('highirect')
% ylabel('lowect')
