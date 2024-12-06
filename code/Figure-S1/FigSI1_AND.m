
clear
clc

load param_ThrePert_SA0SB0_50pc_AND.mat
load param_SameSA0SB0_AND.mat

unp=tri_types(nrr,:);
ind=tri_type_SA0ltSB0;
dir=tri_type_SB0ltSA0;

nr=length(tri_type_SB0ltSA0);

for ii=1:nr

    if tri_type_SA0ltSB0(ii)==1 
        M_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,4);
        H_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,2);
        Tri_Reg_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,1)-SNB_SA0ltSB0(ii,4);
    elseif tri_type_SA0ltSB0(ii)==3
        M_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,4);
        H_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,2);
        Tri_Reg_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,1)-SNB_SA0ltSB0(ii,2);
    elseif tri_type_SA0ltSB0(ii)==2 
        M_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,1)-SNB_SA0ltSB0(ii,4);
        H_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,2);
        Tri_Reg_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,4);
    elseif tri_type_SA0ltSB0(ii)==4 
        M_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,1)-SNB_SA0ltSB0(ii,4);
        H_state_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,2);
        Tri_Reg_SA0ltSB0(ii,1)=SNB_SA0ltSB0(ii,3)-SNB_SA0ltSB0(ii,2);
    end


    if tri_type_SB0ltSA0(ii)==1 
        M_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,4);
        H_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,2);
        Tri_Reg_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,1)-SNB_SB0ltSA0(ii,4);
    elseif tri_type_SB0ltSA0(ii)==3
        M_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,4);
        H_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,2);
        Tri_Reg_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,1)-SNB_SB0ltSA0(ii,2);
    elseif tri_type_SB0ltSA0(ii)==2
        M_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,1)-SNB_SB0ltSA0(ii,4);
        H_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,2);
        Tri_Reg_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,4);
    elseif tri_type_SB0ltSA0(ii)==4
        M_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,1)-SNB_SB0ltSA0(ii,4);
        H_state_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,2);
        Tri_Reg_SB0ltSA0(ii,1)=SNB_SB0ltSA0(ii,3)-SNB_SB0ltSA0(ii,2);
    end

end

ind_types=[sum(ind==1),sum(ind==2),sum(ind==3),sum(ind==4)];
ind_pc=ind_types.*100./sum(ind_types);

dir_types=[sum(dir==1),sum(dir==2),sum(dir==3),sum(dir==4)];
dir_pc=dir_types.*100./sum(dir_types);

unp_types=[sum(unp==1),sum(unp==2),sum(unp==3),sum(unp==4)];
unp_pc=unp_types.*100./sum(unp_types);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the vertices of the triangle
x1 = [0, 300, 300];
y1 = [0, 0, 300];

% Plot the vertices
figure(4)
subplot(3,2,1)
plot(x1, y1, 'k-'); % Plot the lines connecting the vertices

% Fill the area of the triangle
% figure(4)
% subplot(3,3,1)
fill(x1, y1, 'b'); % Fill the triangle with blue color
hold on

% Define the vertices of the triangle
x2 = [0, 0, 300];
y2 = [0, 300, 300];

% Plot the vertices
% figure(4)
% subplot(3,3,1)
plot(x2, y2, 'k-'); % Plot the lines connecting the vertices

% Fill the area of the triangle
% figure(4)
% subplot(3,3,1)
fill(x2, y2, 'g'); % Fill the triangle with blue color
hold on

% figure(4)
% subplot(3,3,1)
plot(SNB_SA0ltSB0(:,1), SNB_SB0ltSA0(:,1))
hold on
plot([0,300],[0,300])

sn1_dir_gt_ind = find(SNB_SB0ltSA0(:,1)>SNB_SA0ltSB0(:,1));
sn1_dir_lt_ind = find(SNB_SB0ltSA0(:,1)<SNB_SA0ltSB0(:,1));
sn1_dir_eq_ind = find(SNB_SB0ltSA0(:,1)==SNB_SA0ltSB0(:,1));
tot_sn1 = length(sn1_dir_gt_ind)+length(sn1_dir_lt_ind)+length(sn1_dir_eq_ind);

sn1_cont=([length(sn1_dir_gt_ind),length(sn1_dir_lt_ind),length(sn1_dir_eq_ind)].*100./tot_sn1)

% Define the vertices of the triangle
x1 = [0, 300, 300];
y1 = [0, 0, 300];

% Plot the vertices
figure(4)
subplot(3,2,3)
plot(x1, y1, 'k-'); % Plot the lines connecting the vertices

% Fill the area of the triangle
% figure(4)
% subplot(2,2,2)
fill(x1, y1, 'b'); % Fill the triangle with blue color
hold on

% Define the vertices of the triangle
x2 = [0, 0, 300];
y2 = [0, 300, 300];

% Plot the vertices
% figure(4)
% subplot(2,2,2)
plot(x2, y2, 'k-'); % Plot the lines connecting the vertices

% Fill the area of the triangle
% figure(4)
% subplot(2,2,2)
fill(x2, y2, 'g'); % Fill the triangle with blue color
hold on

% subplot(2,2,2)
plot(M_state_SA0ltSB0,M_state_SB0ltSA0)
hold on
plot([0,300],[0,300])

m_dir_gt_ind = find(M_state_SB0ltSA0(:,1)>M_state_SA0ltSB0(:,1));
m_dir_lt_ind = find(M_state_SB0ltSA0(:,1)<M_state_SA0ltSB0(:,1));
m_dir_eq_ind = find(M_state_SB0ltSA0(:,1)==M_state_SA0ltSB0(:,1));
tot_m = length(m_dir_gt_ind)+length(m_dir_lt_ind)+length(m_dir_eq_ind);


mstate_cont=([length(m_dir_gt_ind),length(m_dir_lt_ind),length(m_dir_eq_ind)].*100./tot_m)

subplot(3,2,5)
bar([unp_pc;dir_pc;ind_pc]')
% set(gca, 'TickDir', 'out');

text(1:length(unp_pc),unp_pc,num2str(unp_pc'),'vert','bottom','horiz','center','Rotation',90); 

text(1:length(dir_pc),dir_pc,num2str(dir_pc'),'vert','bottom','horiz','center','Rotation',90); 

text(1:length(ind_pc),ind_pc,num2str(ind_pc'),'vert','bottom','horiz','center','Rotation',90); 

set(gca, 'XTickLabel',{'2U2D','1U2D','2U1D','1U1D'})
