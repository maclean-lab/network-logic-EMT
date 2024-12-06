
% This is a code to generate Fig3D

clear
clc

Basin=[];
LocMinVar=[];


% parameters

ksa=23.5891;
ksb=23.5891;
kab=34.6852;
kba=34.6852;
kbb=10.1161;
ab0=478.6441;
ba0=478.6441;
bb0=299.0639;
nsa=6;
nsb=6;
nab=8;
nba=8;
nbb=10;
 ga=0.0491;
 gb=0.0491;
sa0=290.9559;
sb0=290.9559;
ka0=0.4704;
kb0=0.4704;

range = flip([150:10:400]);
SNB = NaN(length(range),4);

for ik = 1:length(range)

    ik

    sa0 = range(ik);

    

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

b1=[0:delb:1000];     % range for argument in potential function
        
        
        S2A = ksa./(1+(s./sa0).^nsa);
        B2A = kba./(1+(b1./ba0).^nba);

        


        a1 = (ka0+S2A+B2A)./ga;

        S2B = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb);
        A2B = kab./(1+(a1./ab0).^nab);
        B2B = kbb.*(b1./bb0).^nbb./(1+(b1./bb0).^nbb);

        
        f = kb0+S2B+A2B+B2B-gb.*b1;

       
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

%SNB(ik,:) = [SN1 SN2 SN3 SN4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Categorizing the response type based on the transition paths

if (SN3>SN1) && (SN1>SN4) && (SN4>SN2)
    disp('2U2D');
    TriType(ik,1) = 22;
    
elseif (SN1>SN3) && (SN3>SN4) && (SN4>SN2)
    disp('1U2D');
    TriType(ik,1) = 12;
    
elseif (SN3>SN1) && (SN1>SN2) && (SN2>SN4)
    disp('2U1D');
    TriType(ik,1) = 21;
   
elseif (SN1>SN3) && (SN3>SN2) && (SN2>SN4)
    disp('1U1D');
    TriType(ik,1) = 11;
    
elseif (SN3>SN4) && (SN4>SN1) && (SN1>SN2)
    disp('DBS');
    TriType(ik,1) = 202;
    
elseif (SN2>SN3)
    SN2 = NaN;
    SN3 = NaN;
    disp('Bistable');
    TriType(ik,1) = 0;
end



SNB(ik,:) = [SN1 SN2 SN3 SN4];

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

end


plot(range,SNB(:,1))
hold on
plot(range,SNB(:,2))
plot(range,SNB(:,3))
plot(range,SNB(:,4))














