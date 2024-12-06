
% This is a code to produce Fig2A, C

clear
clc

Basin=[];
LocMinVar=[];

% parameters

ksa=87.2400;    
ksb=87.2400;    
kab=1;           
kba=1;
kbb=1;
ab0=387.1077;
ba0=387.1077;
bb0=87.2636;
nsa=1;
nsb=1;
nab=3;
nba=3;
nbb=2;
 ga=0.0773;
 gb=0.0773;
sa0=209.7679;
sb0=209.7679;
ka0=0.9486;
kb0=0.9486;



sa0 = 209.7679-0.5*209.7679; % perturbing indirect activation

%sb0 = 209.7679-0.5*209.7679; % perturbing direct activation

%subplot(3,2,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating potential energy landscapes

sig = 0:0.1:300;

nr = length(sig);


BSSS = NaN(nr,3);
BUSS = NaN(nr,3);
ASSS = NaN(nr,3);
AUSS = NaN(nr,3);

for ij=1:nr

s=sig(ij);

sigs(ij,1)=s;   % this is signal (S) which is the bifurcation parameter here



delb=0.1;

b1=[0:delb:3000];     % range for argument in potential function
        
        
        S2A = ksa./(1+(s./sa0).^nsa);
        B2A = kba./(1+(b1./ba0).^nba);

        a1 = (ka0+S2A.*B2A)./ga;

        S2B = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb);
        A2B = kab./(1+(a1./ab0).^nab);
        B2B = kbb.*(b1./bb0).^nbb./(1+(b1./bb0).^nbb);

        
        f = kb0+S2B.*A2B.*B2B-gb.*b1;

       
z=cumtrapz(-f.*delb);         % Potential is the negative of the function integrated

z1=z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the local minima and maxima from the potential landscape

TF_min=islocalmin(z1)';   % finding the local min
LocMinforVar = b1(TF_min==1);
LocMinforPot = z1(TF_min==1);
LocMinforVar2 = a1(TF_min==1);


LocMin_Indx = find(TF_min==1);

TF_max=islocalmax(z1)';   % finding the local max
LocMaxforVar = b1(TF_max==1);
LocMaxforPot = z1(TF_max==1);
LocMaxforVar2 = a1(TF_max==1);

LocMax_Indx = find(TF_max==1);

BSSS(ij,1:length(LocMinforVar)) = LocMinforVar;
BUSS(ij,1:length(LocMaxforVar)) = LocMaxforVar;

ASSS(ij,1:length(LocMinforVar2)) = LocMinforVar2;
AUSS(ij,1:length(LocMaxforVar2)) = LocMaxforVar2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the SSS of B for every S

end

BSSS_sortS=[sigs BSSS(:,1)
            sigs BSSS(:,2)
            sigs BSSS(:,3)];
[~,idx1] = sort(BSSS_sortS(:,2));
        sortedBSSS = BSSS_sortS(idx1,:);

sortedBSSS(any(isnan(sortedBSSS),2),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the USS of B for every S

BUSS_sortS=[sigs BUSS(:,1)
            sigs BUSS(:,2)
            sigs BUSS(:,3)];
[~,idx2] = sort(BUSS_sortS(:,2));
        sortedBUSS = BUSS_sortS(idx2,:);

sortedBUSS(any(isnan(sortedBUSS),2),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the SSS of A for every S

ASSS_sortS=[sigs ASSS(:,1)
            sigs ASSS(:,2)
            sigs ASSS(:,3)];
[~,idx3] = sort(ASSS_sortS(:,2));
        sortedASSS = ASSS_sortS(idx3,:);

sortedASSS(any(isnan(sortedASSS),2),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the USS of A for every S

AUSS_sortS=[sigs AUSS(:,1)
            sigs AUSS(:,2)
            sigs AUSS(:,3)];
[~,idx4] = sort(AUSS_sortS(:,2));
        sortedAUSS = AUSS_sortS(idx4,:);

sortedAUSS(any(isnan(sortedAUSS),2),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the SN points from the USSS (large gaps in the data for USSS of B)

diffBUSS2=diff(sortedBUSS(:,2));
maxBUSS2 = max(diffBUSS2);
IndxUSS2=find(diffBUSS2==maxBUSS2);

SN_B=[sortedBUSS(1,1) sortedBUSS(1,2);
      sortedBUSS(IndxUSS2,1) sortedBUSS(IndxUSS2,2);
      sortedBUSS(IndxUSS2+1,1) sortedBUSS(IndxUSS2+1,2);
      sortedBUSS(end,1) sortedBUSS(end,2)];

SN1=SN_B(1);
SN2=SN_B(2);
SN3=SN_B(3);
SN4=SN_B(4);

%SNB(1,:) = [SN1 SN2 SN3 SN4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Categorizing the response type based on the transition paths

if (SN3>SN1) && (SN1>SN4) && (SN4>SN2)
    disp('2U2D')
    
    
elseif (SN1>SN3) && (SN3>SN4) && (SN4>SN2)
    disp('1U2D')
    
    
elseif (SN3>SN1) && (SN1>SN2) && (SN2>SN4)
    disp('2U1D')
   
   
elseif (SN1>SN3) && (SN3>SN2) && (SN2>SN4)
    disp('1U1D')
    
    
elseif (SN3>SN4) && (SN4>SN1) && (SN1>SN2)
    disp('DBS')
    

elseif (SN2>SN3)
    SN3 = NaN;
    SN4 = NaN;
    disp('Bistable');
    
end



diffAUSS2=diff(sortedAUSS(:,2));
maxAUSS2 = max(diffAUSS2);
IndxAUSS2=find(diffAUSS2==maxAUSS2);

SN_A=[sortedAUSS(1,1) sortedAUSS(1,2);
      sortedAUSS(IndxAUSS2,1) sortedAUSS(IndxAUSS2,2);
      sortedAUSS(IndxAUSS2+1,1) sortedAUSS(IndxAUSS2+1,2);
      sortedAUSS(end,1) sortedAUSS(end,2)];

% separating the lower and upper unstable SS
diff_data = diff(sortedBUSS(:,2));
[maxValue, index] = max(diff_data);
sortedBUSS_Low = sortedBUSS(1:index,:);
sortedBUSS_Up = sortedBUSS(index+1:end,:);

%subplot(2,2,1)
plot(sortedBSSS(:,1),sortedBSSS(:,2),'black.')
hold on
plot(sortedBUSS_Low(:,1),sortedBUSS_Low(:,2),'b--')
hold on
plot(sortedBUSS_Up(:,1),sortedBUSS_Up(:,2),'b--')
hold on
%plot(sortedBUSS(:,1),sortedBUSS(:,2),'green.')
% hold on
plot(SN_B(1,1),SN_B(1,2),'o')
plot(SN_B(2,1),SN_B(2,2),'o')
plot(SN_B(3,1),SN_B(3,2),'o')
plot(SN_B(4,1),SN_B(4,2),'o')
hold on

% Set axis labels
xlabel('S');
ylabel('B');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To plot the loci of the SN points as a function of sa0 (indirect
% activation strength) or sb0 (direct activation strength), rerun this code
% for a range of either sa0 or sb0 and record the SN points at each sa0 or
% sb0. The loci of the SN points at different SA0/SB0 is stored in the
% matlab file 'SNpoints_VaryingSA0_AND.mat' and 'SNpoints_VaryingSB0_AND.mat'
% Call these files to the matlab code 'SNpoints_VaryingSA0SB0_plots_AND.m' to get
% figures 2B,D. 














