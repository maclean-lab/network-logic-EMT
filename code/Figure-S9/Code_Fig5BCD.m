clear
clc

% This MATLAB code is to plot figures Fig 5B, C, D

%load stochEM_sigma1_Unpert_OR.mat
%load stochEM_sigma1_25pcInd_OR.mat
%load stochEM_sigma1_25pcDir_OR.mat

%load stochEM_sigma1_Unpert_AND.mat
%load stochEM_sigma1_25pcInd_AND.mat
%load stochEM_sigma1_25pcDir_AND.mat

%load stochEM_sigma1_V1_Unpert_AND.mat
%load stochEM_sigma1_25pcInd_AND.mat
%load stochEM_sigma1_25pcDir_AND.mat

res_aa = A_sig;

res_bb = B_sig;

max_aa = max(res_aa);
min_aa = min(res_aa);
max_bb = max(res_bb);
min_bb = min(res_bb);


for ij=1:8
    ii=ij;

    nr=1000;


res_aa1 = (res_aa-min_aa)./(max_bb-min_aa);
res_bb1 = (res_bb-min_bb)./(max_bb-min_bb);


x=(res_aa1((ii-1)*nr+1:ii*nr));
y=(res_bb1((ii-1)*nr+1:ii*nr));


data_aa(1:nr,ii)=x;
data_bb(1:nr,ii)=y;


% Estimate density for each data point
densityValues = ksdensity([x y], [x y]);

% Ensure densityValues is a column vector with the same length as x and y
densityValues = densityValues(:);

% % Create scatter plot with color code for data density
figure(8);
%subplot(3,8,ii)
subplot(4,4,ii)
%subplot(5,5,ii+10)
scatter(x, y, 10, densityValues, 'filled');
colormap('jet'); 
%colorbar;
hold on


xrange=[0,1];
yrange=[0,1];


% Set axis labels and title
xlabel('A');
ylabel('B');
xlim(xrange)
ylim(yrange)

end

