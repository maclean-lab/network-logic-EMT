
clear
clc

load varSDE_Multi_L3_N0o1_Fp0o066_MT.mat

time_points = [5,15,30,50];

ind = find(ismember(t, time_points));

%var(1) is Gata6; var(2) is Nanog

nr = 20;


%figure('Position', [100, 100, 1000, 800]);

ICM_G1=NaN(25,20);
ICM_N1=NaN(25,20);
Plike_G1=NaN(25,20);
Plike_N1=NaN(25,20);
Elike_G1=NaN(25,20);
Elike_N1=NaN(25,20);
PrE_G1=NaN(25,20);
PrE_N1=NaN(25,20);
Epi_G1=NaN(25,20);
Epi_N1=NaN(25,20);

ICM_G2=NaN(25,20);
ICM_N2=NaN(25,20);
Plike_G2=NaN(25,20);
Plike_N2=NaN(25,20);
Elike_G2=NaN(25,20);
Elike_N2=NaN(25,20);
PrE_G2=NaN(25,20);
PrE_N2=NaN(25,20);
Epi_G2=NaN(25,20);
Epi_N2=NaN(25,20);

ICM_G3=NaN(25,20);
ICM_N3=NaN(25,20);
Plike_G3=NaN(25,20);
Plike_N3=NaN(25,20);
Elike_G3=NaN(25,20);
Elike_N3=NaN(25,20);
PrE_G3=NaN(25,20);
PrE_N3=NaN(25,20);
Epi_G3=NaN(25,20);
Epi_N3=NaN(25,20);

ICM_G4=NaN(25,20);
ICM_N4=NaN(25,20);
Plike_G4=NaN(25,20);
Plike_N4=NaN(25,20);
Elike_G4=NaN(25,20);
Elike_N4=NaN(25,20);
PrE_G4=NaN(25,20);
PrE_N4=NaN(25,20);
Epi_G4=NaN(25,20);
Epi_N4=NaN(25,20);

for ii=1:nr

    var = cell2mat(var_values(ii,1));

    nrr = size(var,3);

    cnt1(ii,1) = 0;     % ICM
    cnt1(ii,2) = 0;     % PrE like
    cnt1(ii,3) = 0;     % Epi like
    cnt1(ii,4) = 0;     % PrE
    cnt1(ii,5) = 0;     % Epi

    cnt2(ii,1) = 0;     % ICM
    cnt2(ii,2) = 0;     % PrE like
    cnt2(ii,3) = 0;     % Epi like
    cnt2(ii,4) = 0;     % PrE
    cnt2(ii,5) = 0;     % Epi

    cnt3(ii,1) = 0;     % ICM
    cnt3(ii,2) = 0;     % PrE like
    cnt3(ii,3) = 0;     % Epi like
    cnt3(ii,4) = 0;     % PrE
    cnt3(ii,5) = 0;     % Epi

    cnt4(ii,1) = 0;     % ICM (black)
    cnt4(ii,2) = 0;     % PrE like  (cyan)
    cnt4(ii,3) = 0;     % Epi like  (magenta)
    cnt4(ii,4) = 0;     % PrE   (blue)
    cnt4(ii,5) = 0;     % Epi   (red)
    

    for ij=1:nrr

       
        
        tp = ind(1); % time point 1
        if var(tp,1,ij) <= 0.5 && var(tp,2,ij) <= 0.5
            ICM_G1(ij,ii) = var(tp,1,ij);
            ICM_N1(ij,ii) = var(tp,2,ij);
            cnt1(ii,1) = cnt1(ii,1)+1;
        elseif (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1) && (var(tp,2,ij) <= 0.5)
            Plike_G1(ij,ii) = var(tp,1,ij);
            Plike_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,2) = cnt1(ii,2)+1;
        elseif (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1) && (var(tp,1,ij) <= 0.5)
            Elike_G1(ij,ii) = var(tp,1,ij);
            Elike_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,3) = cnt1(ii,3)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) <= 0.5)
            PrE_G1(ij,ii) = var(tp,1,ij);
            PrE_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,4) = cnt1(ii,4)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) <= 0.5)
            Epi_G1(ij,ii) = var(tp,1,ij);
            Epi_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,5) = cnt1(ii,5)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1)
            Plike_G1(ij,ii) = var(tp,1,ij);
            Plike_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,2) = cnt1(ii,2)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1)
            Elike_G1(ij,ii) = var(tp,1,ij);
            Elike_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,3) = cnt1(ii,3)+1;
        else 
            ICM_G1(ij,ii) = var(tp,1,ij);
            ICM_N1(ij,ii) = var(tp,2,ij);
            
            cnt1(ii,1) = cnt1(ii,1)+1;
        end

        xlim([0 2])
        ylim([0 2])
        xlabel('Gata6')
        ylabel('Nanog')
        title('TimePoint = 5')

        tp = ind(2); % time point 2
        if var(tp,1,ij) <= 0.5 && var(tp,2,ij) <= 0.5
            ICM_G2(ij,ii) = var(tp,1,ij);
            ICM_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,1) = cnt2(ii,1)+1;
        elseif (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1) && (var(tp,2,ij) <= 0.5)
            Plike_G2(ij,ii) = var(tp,1,ij);
            Plike_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,2) = cnt2(ii,2)+1;
        elseif (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1) && (var(tp,1,ij) <= 0.5)
            Elike_G2(ij,ii) = var(tp,1,ij);
            Elike_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,3) = cnt2(ii,3)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) <= 0.5)
            PrE_G2(ij,ii) = var(tp,1,ij);
            PrE_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,4) = cnt2(ii,4)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) <= 0.5)
            Epi_G2(ij,ii) = var(tp,1,ij);
            Epi_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,5) = cnt2(ii,5)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1)
            Plike_G2(ij,ii) = var(tp,1,ij);
            Plike_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,2) = cnt2(ii,2)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1)
            Elike_G2(ij,ii) = var(tp,1,ij);
            Elike_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,3) = cnt2(ii,3)+1;
        else 
            ICM_G2(ij,ii) = var(tp,1,ij);
            ICM_N2(ij,ii) = var(tp,2,ij);
            
            cnt2(ii,1) = cnt2(ii,1)+1;
        end
        xlim([0 2])
        ylim([0 2])
        xlabel('Gata6')
        ylabel('Nanog')
        title('TimePoint = 15')


        tp = ind(3); % time point 2
        if var(tp,1,ij) <= 0.5 && var(tp,2,ij) <= 0.5
            ICM_G3(ij,ii) = var(tp,1,ij);
            ICM_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,1) = cnt3(ii,1)+1;
        elseif (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1) && (var(tp,2,ij) <= 0.5)
            Plike_G3(ij,ii) = var(tp,1,ij);
            Plike_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,2) = cnt3(ii,2)+1;
        elseif (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1) && (var(tp,1,ij) <= 0.5)
            Elike_G3(ij,ii) = var(tp,1,ij);
            Elike_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,3) = cnt3(ii,3)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) <= 0.5)
            PrE_G3(ij,ii) = var(tp,1,ij);
            PrE_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,4) = cnt3(ii,4)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) <= 0.5)
            Epi_G3(ij,ii) = var(tp,1,ij);
            Epi_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,5) = cnt3(ii,5)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1)
            Plike_G3(ij,ii) = var(tp,1,ij);
            Plike_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,2) = cnt3(ii,2)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1)
            Elike_G3(ij,ii) = var(tp,1,ij);
            Elike_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,3) = cnt3(ii,3)+1;
        else 
            ICM_G3(ij,ii) = var(tp,1,ij);
            ICM_N3(ij,ii) = var(tp,2,ij);
            
            cnt3(ii,1) = cnt3(ii,1)+1;
        end
        xlim([0 2])
        ylim([0 2])
        xlabel('Gata6')
        ylabel('Nanog')
        title('TimePoint = 30')



        tp = ind(4); % time point 2
        if var(tp,1,ij) <= 0.5 && var(tp,2,ij) <= 0.5
            ICM_G4(ij,ii) = var(tp,1,ij);
            ICM_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,1) = cnt4(ii,1)+1;
        elseif (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1) && (var(tp,2,ij) <= 0.5)
            Plike_G4(ij,ii) = var(tp,1,ij);
            Plike_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,2) = cnt4(ii,2)+1;
        elseif (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1) && (var(tp,1,ij) <= 0.5)
            Elike_G4(ij,ii) = var(tp,1,ij);
            Elike_N4(ij,ii) = var(tp,2,ij);
           
            cnt4(ii,3) = cnt4(ii,3)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) <= 0.5)
            PrE_G4(ij,ii) = var(tp,1,ij);
            PrE_N4(ij,ii) = var(tp,2,ij);
           
            cnt4(ii,4) = cnt4(ii,4)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) <= 0.5)
            Epi_G4(ij,ii) = var(tp,1,ij);
            Epi_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,5) = cnt4(ii,5)+1;
        elseif (var(tp,1,ij) > 1) && (var(tp,2,ij) > 0.5 && var(tp,2,ij) <= 1)
            Plike_G4(ij,ii) = var(tp,1,ij);
            Plike_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,2) = cnt4(ii,2)+1;
        elseif (var(tp,2,ij) > 1) && (var(tp,1,ij) > 0.5 && var(tp,1,ij) <= 1)
            Elike_G4(ij,ii) = var(tp,1,ij);
            Elike_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,3) = cnt4(ii,3)+1;
        else 
            ICM_G4(ij,ii) = var(tp,1,ij);
            ICM_N4(ij,ii) = var(tp,2,ij);
            
            cnt4(ii,1) = cnt4(ii,1)+1;
        end
        xlim([0 2])
        ylim([0 2])
        xlabel('Gata6')
        ylabel('Nanog')
        title('TimePoint = 50')


        %sgtitle(sprintf('SDE (Multiplicative) for Logic 3; Noise Strength = 0.1; MT; Fp = 0o066'));
    end

end

subplot(8,9,24)
scatter(ICM_G1(:),ICM_N1(:),'black','filled')
hold on
scatter(Plike_G1(:),Plike_N1(:),'cyan','filled')
scatter(Elike_G1(:),Elike_N1(:),'magenta','filled')
scatter(PrE_G1(:),PrE_N1(:),'blue','filled')
scatter(Epi_G1(:),Epi_N1(:),'red','filled')
xlim([0 2])
ylim([0 2])
xlabel('Gata6')
ylabel('Nanog')
title('TimePoint = 5')


subplot(8,9,25)
scatter(ICM_G2(:),ICM_N2(:),'black','filled')
hold on
scatter(Plike_G2(:),Plike_N2(:),'cyan','filled')
scatter(Elike_G2(:),Elike_N2(:),'magenta','filled')
scatter(PrE_G2(:),PrE_N2(:),'blue','filled')
scatter(Epi_G2(:),Epi_N2(:),'red','filled')
xlim([0 2])
ylim([0 2])
xlabel('Gata6')
ylabel('Nanog')
title('TimePoint = 15')

subplot(8,9,26)
scatter(ICM_G3(:),ICM_N3(:),'black','filled')
hold on
scatter(Plike_G3(:),Plike_N3(:),'cyan','filled')
scatter(Elike_G3(:),Elike_N3(:),'magenta','filled')
scatter(PrE_G3(:),PrE_N3(:),'blue','filled')
scatter(Epi_G3(:),Epi_N3(:),'red','filled')
xlim([0 2])
ylim([0 2])
xlabel('Gata6')
ylabel('Nanog')
title('TimePoint = 30')

subplot(8,9,27)
scatter(ICM_G4(:),ICM_N4(:),'black','filled')
hold on
scatter(Plike_G4(:),Plike_N4(:),'cyan','filled')
scatter(Elike_G4(:),Elike_N4(:),'magenta','filled')
scatter(PrE_G4(:),PrE_N4(:),'blue','filled')
scatter(Epi_G4(:),Epi_N4(:),'red','filled')
xlim([0 2])
ylim([0 2])
xlabel('Gata6')
ylabel('Nanog')
title('TimePoint = 50')



cnt11 = sum(cnt1);
cnt22 = sum(cnt2);
cnt33 = sum(cnt3);
cnt44 = sum(cnt4);

colors = [
    0 0 0;    % Black
    0 1 1;    % Cyan
    1 0 1;    % Magenta
    0 0 1;    % Blue
    1 0 0     % Red
];



% subplot(4,4,5)
% bar(1:length(cnt11), cnt11*100/sum(cnt11), 'FaceColor', 'flat', 'CData', colors);
% 
% xticks(1:length(cnt11)); % Set tick positions at 1, 2, 3, 4, 5
% xticklabels({'ICM', 'PrE like', 'Epi like', 'PrE', 'Epi'}); % Set custom labels
% 
% ylim([0 100])
% 
% ylabel('Percent of cells')
% title('TimePoint = 5')
% 
% 
% subplot(4,4,6)
% bar(1:length(cnt22), cnt22*100/sum(cnt22), 'FaceColor', 'flat', 'CData', colors);
% 
% xticks(1:length(cnt22)); % Set tick positions at 1, 2, 3, 4, 5
% xticklabels({'ICM', 'PrE like', 'Epi like', 'PrE', 'Epi'}); % Set custom labels
% 
% ylim([0 100])
% 
% ylabel('Percent of cells')
% title('TimePoint = 15')
% 
% subplot(4,4,7)
% bar(1:length(cnt33), cnt33*100/sum(cnt33), 'FaceColor', 'flat', 'CData', colors);
% 
% xticks(1:length(cnt33)); % Set tick positions at 1, 2, 3, 4, 5
% xticklabels({'ICM', 'PrE like', 'Epi like', 'PrE', 'Epi'}); % Set custom labels
% 
% ylim([0 100])
% 
% ylabel('Percent of cells')
% title('TimePoint = 30')
% 
% subplot(4,4,8)
% bar(1:length(cnt44), cnt44*100/sum(cnt44), 'FaceColor', 'flat', 'CData', colors);
% 
% xticks(1:length(cnt44)); % Set tick positions at 1, 2, 3, 4, 5
% xticklabels({'ICM', 'PrE like', 'Epi like', 'PrE', 'Epi'}); % Set custom labels
% 
% ylim([0 100])
% 
% ylabel('Percent of cells')
% title('TimePoint = 50')

subplot(8,9,33)
pie(cnt11)
colormap(colors);

subplot(8,9,34)
pie(cnt22)
colormap(colors);

subplot(8,9,35)
pie(cnt33)
colormap(colors);

subplot(8,9,36)
pie(cnt44)
colormap(colors);

