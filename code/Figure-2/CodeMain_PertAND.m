
clear
clc

% load the tristable parameter sets generated from the
% CodeMain_RandParSearch_AND.m
% After loading it here, systematically perturb the direct/indirect
% regulation by x% and rerun the code to see how many of the original
% parameter sets generates tristability
% carefully save these tristable parameter sets obtained after the
% perturbation
% These parameter sets are stored in param_ThrePert_SA0SB0_25pc_AND.mat

tic

% counter for bifurcation types
Type1=0;    
Type2=0;
Type3=0;
Type4=0;
Type5=0;
Type6=0;
Type7=0;

load param_SameSA0SB0_AND.mat
param = param_SameSA0SB0;

nr = length(param(1,:));

for kk=1:nr                 % loop for parameters
    kk;

ksa=param(1,kk);
ksb=param(2,kk);
kab=param(3,kk);
kba=param(4,kk);
kbb=param(5,kk);
ab0=param(6,kk);
ba0=param(7,kk);
bb0=param(8,kk);
nsa=param(9,kk);
nsb=param(10,kk);
nab=param(11,kk);
nba=param(12,kk);
nbb=param(13,kk);
 ga=param(14,kk);
 gb=param(15,kk);
sa0=param(16,kk);
sb0=param(17,kk);
ka0=param(18,kk);
kb0=param(19,kk);

%sa0 = param(16,kk)-0.25*param(16,kk);   % Indirect perturbation

%sb0 = param(17,kk)-0.25*param(17,kk);   % Direct perturbation


NN1=300;
NN2=3000;
del1=01;
del2=0.1;
int1=0.1;
int2=0.1;
intf=2500;


nmbr_peaks_ss=[];
nmbr_peaks_us=[];
tri_par=[];
bfd=[];
yx=[];

ic1=0;


    % *********************************************************************
    % Finding no. of stable and unstable steady states with variation of
    % bifurcation parameter. Stable and unstable steady states are given 
    % by minima and maxima in the PEF respectively. 
    % *********************************************************************
    for jj=1:NN1             % loop for bifurcation parameter
        ik(jj)=jj.*del1;
        jj;

       s=ik(jj);

        b1=[0:int1:intf];     % range for argument in potential function

        
        
        S2A = ksa./(1+(s./sa0).^nsa);
        B2A = kba./(1+(b1./ba0).^nba);

        a1 = (ka0+S2A.*B2A)./ga;

        S2B = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb);
        A2B = kab./(1+(a1./ab0).^nab);
        B2B = kbb.*(b1./bb0).^nbb./(1+(b1./bb0).^nbb);
        
        
        f = kb0+S2B.*A2B.*B2B-gb.*b1;

        z=cumtrapz(-f);         % z: potential function and it is 
                                % the negative of the function (f) integrated
        z1=z-min(z);            % rescaling of potential against the global
                                % minima of the pontential function
        z2=-z1;
        % figure(2)
        % plot(x,z2)
        % findpeaks(z2)
        nmbr_peaks_ss(jj)=numel(findpeaks(z2));     % no. of minima in PEF
        nmbr_peaks_us(jj)=numel(findpeaks(z1));     % no. of mamima in PEF

    end
    % *********************************************************************
    % End of finding no. of steady states calculation
    % *********************************************************************
    
    % *********************************************************************
    % Calculation of bifurcation diagram for the parameter combination that
    % generates more than one stable steady states(only for multistability)
    % *********************************************************************
    
    % Counting the no. of stable steady state
    if (nmbr_peaks_ss(1)>1)         % avoids irreversible switch
        continue
    elseif (nmbr_peaks_us(end)>=1)
        continue
    elseif (nmbr_peaks_us(1)>=1)
        continue
    elseif (nmbr_peaks_ss(end)>1)   % restricts bifurcation analysis beyond
         continue                   % the bifurcation parameter range
    elseif any(nmbr_peaks_ss==3)
    %    disp('Tristability found')
        ic1=ic1+1;
        tri_par(ic1)=kk;            % index for the parameter combination
        kk                         % that results tristability

        bfd=ones(NN2,6)*NaN;

        for j=1:NN2      % Loop for bifurcation parameter
            n(j)=j*del2;
            j;

            s=n(j);

            
            
           b2=[0:int1:intf];
            
        
            S2A = ksa./(1+(s./sa0).^nsa);
            B2A = kba./(1+(b2./ba0).^nba);

            a2 = (ka0+S2A.*B2A)./ga;

            S2B = ksb.*(s./sb0).^nsb./(1+(s./sb0).^nsb);
            A2B = kab./(1+(a2./ab0).^nab);
            B2B = kbb.*(b2./bb0).^nbb./(1+(b2./bb0).^nbb);
            
            
            f = kb0+S2B.*A2B.*B2B-gb.*b2;
               
            z=cumtrapz(-f);
            z1=z-min(z);
            z2=z1*-1;      

            %                     findpeaks(z2)

            % nmbr(j)=numel(findpeaks(z2));  % gives total number of peaks

            [pks_stb locs_stb] = findpeaks(z2); % pks_stb: no. of stable ss
            rz_stb=repmat(n(j),size(locs_stb)); % locs_stb: location of ss
            dz_stb=[rz_stb;locs_stb.*int2]';     % storing the coordinates of 
                                                % stable ss.

            % rearranging the data for stable manifold 
            bfd(j,1)=n(j);
            if(size(dz_stb,1)==3)
                bfd(j,2)=dz_stb(1,2);
                bfd(j,3)=NaN;
                bfd(j,4)=dz_stb(2,2);
                bfd(j,5)=NaN;
                bfd(j,6)=dz_stb(3,2);
            elseif(size(dz_stb,1)==2)
                bfd(j,2)=dz_stb(1,2);
                bfd(j,3)=NaN;
                bfd(j,4)=dz_stb(2,2);
                bfd(j,5)=NaN;
                bfd(j,6)=NaN;
            elseif(size(dz_stb,1)==1)
                bfd(j,2)=dz_stb(1,2);
                bfd(j,3)=NaN;
                bfd(j,4)=NaN;
                bfd(j,5)=NaN;
                bfd(j,6)=NaN;
            else
                bfd(j,2)=NaN;
                bfd(j,3)=NaN;
                bfd(j,4)=NaN;
                bfd(j,5)=NaN;
                bfd(j,6)=NaN;
            end

            [pks_us locs_us] = findpeaks(z1);   % pks_stb: no. of unstable ss
            rz_us=repmat(n(j),size(locs_us));   % locs_stb: location of ss
            dz_us=[rz_us;locs_us.*int2]';       % storing the coordinates of 
                                                % unstable ss.
            if(size(pks_us,2)==2)
                bfd(j,3)=dz_us(1,2);
                bfd(j,5)=dz_us(2,2);
            elseif (size(pks_us,2)==1)
                bfd(j,3)=dz_us(1,2);
                bfd(j,5)=NaN;
            end
        end

        yx=[bfd(:,1) bfd(:,2)
            bfd(:,1) bfd(:,3)
            bfd(:,1) bfd(:,4)
            bfd(:,1) bfd(:,5)
            bfd(:,1) bfd(:,6)];
        [~,idx1] = sort(yx(:,2));
        bdia = yx(idx1,:);
        
        bfr=bdia(bdia(:,2)>0);
        steady=(bdia(1:size(bfr,1),:));

%         plot of bifurcation diagram (only for multistability)
        figure(12)
        plot(bdia(:,1),bdia(:,2))
        hold on

        % finding the saddle-node bifurcation point
        unsta=[bfd(:,1) bfd(:,3)
            bfd(:,1) bfd(:,5)];
        [~,idx2] = sort(unsta(:,2));
        yy1 = unsta(idx2,:);
                        figure(12)
                        plot(yy1(:,1),yy1(:,2),'.')
                        hold on
        yz=yy1(yy1(:,2)>0);
        yz1=(yy1(1:size(yz,1),:));

        %                 if (size(yz,1)==0)
        %                     disp('monostable')
        %                     Type1=Type1+1;
        %                     continue
        %                 end

        x1=yz(1);
        y1=yz1(1,2);
        x4=yz(end);
        y4=yz1(end,2);

        diffr=diff(yz1(:,2));

        [~,ind]=max(diffr);
        x2=yz1(ind,1);
        y2=yz1(ind,2);

        x3=yz1(ind+1,1);
        y3=yz1(ind+1,2);
        
        yx1=yx(:,2);
        yx2=yx1(yx1>0);
        
        maxst=max(steady(:,1));
        index_maxst=find(steady(:,1)==maxst);
        yvalmax=steady(index_maxst,2);
        
        minst=min(steady(:,1));
        index_minst=find(steady(:,1)==minst);
        yvalmin=steady(index_minst,2);
        
        maxy=max(steady(:,2));
        index_maxy=find(steady(:,2)==maxy);
        xvalmaxy=steady(index_maxy(1),1);
        
        miny=min(steady(:,2));
        index_miny=find(steady(:,2)==miny);
        xvalminy=steady(index_miny(1),1);
        
        
        if (xvalmaxy<xvalminy)
           sn1=x4;
           ybp1=y4;
           sn4=x1;
           ybp4=y1;
           sn3=x2;
           ybp3=y2;
           sn2=x3;
           ybp2=y3;
        elseif (xvalmaxy>xvalminy)
           sn4=x4;
           ybp4=y4;
           sn1=x1;
           ybp1=y1;
           sn2=x2;
           ybp2=y2;
           sn3=x3;
           ybp3=y3;
        end

        if (sn2==sn4)
            sn4=sn4+(-del2+(del2+del2).*rand(1));
        end

        if (sn1==sn3)
            sn3=sn3+(-del2+(del2+del2).*rand(1));
        end

        if (sn2>sn3)
            sn2=NaN;
            sn3=NaN;
            y2=NaN;
            y3=NaN;
	    Type2=Type2+1;
	    disp('BS')
        end

        bp1=[sn1 ybp1];
        bp2=[sn2 ybp2];
        bp3=[sn3 ybp3];
        bp4=[sn4 ybp4];
        bifr_pts=[bp1; bp2; bp3; bp4];

        SNB(ic1,:)=[sn1,sn2,sn3,sn4];

        % determination of type of bifurcation in tistability
        if (sn1>sn3) && (sn3>sn2) && (sn2>sn4)
            
            disp('1U1D')
            tri_types(ic1,1)=4;
            
        elseif (sn3>sn1) && (sn1>sn2) && (sn2>sn4)
            
            disp('2U1D')
            tri_types(ic1,1)=3;
            
        elseif (sn1>sn3) && (sn3>sn4) && (sn4>sn2)
            
            disp('1U2D')
            tri_types(ic1,1)=2;
            
        elseif (sn3>sn1) && (sn1>sn4) && (sn1>sn2)
            
            disp('2U2D')
            tri_types(ic1,1)=1;
           
        elseif (sn3>sn4) && (sn4>sn1) && (sn1>sn2)
            
            disp('DBS')
            tri_types(ic1,1)=5;
        elseif (sn1==sn2)
            
            disp('BS')
            tri_types(ic1,1)=0;
        elseif (sn3==sn4)
           
            disp('BS')
            tri_types(ic1,1)=0;
        end

%         plot(bp1(:,1),bp1(:,2),'r*')
%         plot(bp2(:,1),bp2(:,2),'b*')
%         plot(bp3(:,1),bp3(:,2),'g*')
%         plot(bp4(:,1),bp4(:,2),'magenta*')
%         hold on
                
    end

end

           
sss=[bdia(:,1) bdia(:,2)];
uss=[yy1(:,1) yy1(:,2)];


toc        
        
