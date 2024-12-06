
% MATLAB code to generate bifurcation diagrams using effective potential
% energy method


clear
clc

Basin=[];
LocMinVar=[];

for ik = 1:6

    ik
    
% Parameter input
ka0=0;          % kA0 in manuscript
kb0=0;          % KB0 in manuscript
ksa=8;          % kSA in manuscript
ksb=8;          % kSB in manuscript
kab=10;         % kAB in manuscript
kba=10;         % kBA in manuscript
kbb=7;          % kBB in manuscript
sa0=230;        % SA0 in manuscript
sb0=230;        % SB0 in manuscript
ab0=80;         % AB0 in manuscript
ba0=80;         % BA0 in manuscript
bb0=36;         % BB0 in manuscript
nsa=4;          % nSA in manuscript
nsb=4;          % nSB in manuscript
nab=6;          % nAB in manuscript
nba=6;          % nBA in manuscript
nbb=3;          % nBB in manuscript
ga=0.13;        % gammaA in manuscript
gb=0.13;        % gammaB in manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change these parameters to generate different response type


    if ik==1
        ba0 = 88;   %2U2D
    elseif ik==2
        ba0 = 80;   %2U1D
    elseif ik==3
        nbb = 5;    %1U2D
    elseif ik==4
        nbb = 4;    %1U1D
    elseif ik==5
        ba0 = 95;   %DBS
    elseif ik==6
        bb0 = 50;     %BS
    end


subplot(2,3,ik)

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

b1=[0:delb:600];     % range for argument in potential function
        
        
        S2A = ksa./(1+(s./sa0).^nsa);
        B2A = kba./(1+(b1./ba0).^nba);

        a1 = (ka0+S2A+B2A)./ga;

        S2B = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb);
        A2B = kab./(1+(a1./ab0).^nab);
        B2B = kbb.*(b1./bb0).^nbb./(1+(b1./bb0).^nbb);

        
        f = kb0+S2B+A2B+B2B-gb.*b1;

       
z=cumtrapz(-f.*delb); % Potential is the negative of the function integrated
z1=z-min(z);

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


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the SSS of B for every S and sorting them

BSSS_sortS=[sigs BSSS(:,1)
            sigs BSSS(:,2)
            sigs BSSS(:,3)];
[~,idx1] = sort(BSSS_sortS(:,2));
        sortedBSSS = BSSS_sortS(idx1,:);

sortedBSSS(any(isnan(sortedBSSS),2),:) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the USS of B for every S and sorting them

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
    SN_B(2,1) = NaN;
    SN_B(3,1) = NaN;
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
fnd = find(diff(sortedBUSS(:,2))>1);
sortedBUSS_Low = sortedBUSS(1:fnd,:);
sortedBUSS_Up = sortedBUSS(fnd+1:end,:);

%subplot(2,2,1)
plot(sortedBSSS(:,1),sortedBSSS(:,2),'black.')
hold on
plot(sortedBUSS_Low(:,1),sortedBUSS_Low(:,2),'black--')
hold on
plot(sortedBUSS_Up(:,1),sortedBUSS_Up(:,2),'black--')
hold on

if ik==6
plot(sortedBUSS(:,1),sortedBUSS(:,2),'black--')
hold on
end

plot(SN_B(1,1),SN_B(1,2),'o','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880])
plot(SN_B(2,1),SN_B(2,2),'o','MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])
plot(SN_B(3,1),SN_B(3,2),'o','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1])
plot(SN_B(4,1),SN_B(4,2),'o','MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[1 0 1])


% Set axis labels
xlabel('S');
ylabel('B');

end
















