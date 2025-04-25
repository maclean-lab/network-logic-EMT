
clear
clc

load param_ThrePert_SA0SB0_25pc_OR.mat
%load param_SameSA0SB0_AND.mat

EMT_Win_Ind = SNB_SA0ltSB0(:,1) - SNB_SA0ltSB0(:,4);
EMT_Win_Dir = SNB_SB0ltSA0(:,1) - SNB_SB0ltSA0(:,4);

figure(3)
subplot(2,2,4)
plot(EMT_Win_Ind,EMT_Win_Dir,'o')
hold on
plot([0,300],[0,300])

DirGtInd = length(find(EMT_Win_Dir > EMT_Win_Ind))*100/length(EMT_Win_Ind)

IndrGtDir = 100-DirGtInd