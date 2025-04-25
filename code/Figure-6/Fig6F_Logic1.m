
% This code is to generate Figures 6F



  clear
  clc

    tspan = [0 20];
    c0 = [0 0 2.8 0.1 0.0619]';  % G, N, FR, ERK (0.25), Fs; TriHalf = 0.0619
    [t,c] = ode45(@(t,c) Logic1(t,c), tspan, c0);

    figure(1)
    subplot(2,2,2)
    plot(t,c(:,1),'blue-')
    hold on
    plot(t,c(:,2),'red-')
    ylim([0 2])
    %title('Logic 1, WT, Fp = 0.0619 (TriHalf)');

function dcdt = Logic1(t,c)

  dcdt = zeros(5,1);

  G1 = c(1);
  N1 = c(2);
  FR1 = c(3);
  ERK1 = c(4);
  Fs1 = c(5);
 
 
  % Parameters
Kig = 2; vsg1 = 1.202; vsg2 = 1; Kag1 = 0.28; Kag2 = 0.55; kdg = 1;
Kin2 = 2; vsn1 = 0.856; Kin1 = 0.28; vsn2 = 1; Kan = 0.55; kdn = 1;
vsfr1 = 2.8; Kifr = 0.5; vsfr2 = 2.8; Kafr = 0.5; kdfr = 1;
va = 20; Kd = 2; Ka = 0.7; vi = 3.3; Ki = 0.7;
vex = 0; vsf = 0.6; Kaf = 5; kdf = 0.09;
q = 4; r = 3; s = 4; w = 4; u = 3; v = 4; z = 4;


 Fp1 = Fs1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % for gata6 mutant
 %vsg1 = 0; vsg2 = 0;

 % Phase 1: Nanog is not expressed when fgf4 is administered early
 % vex = 1 at t=0;

  %vex = 1;

 % Phase 2: Nanog is maintained when fgf4 is given at later time point
 % vex = 1 at t=3;

 % if t>=10
 % 
 %    vex = 1;
 % 
 % else
 %     vex = 0;
 % end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % Nanog mutant
 vsn1 = 0; vsn2 = 0;

 %cell1: no inhibitor (va=20)
 % cell 2: inhibitor given at t=1 (Phase 1)
 % cell 3: inhibitor given at t = 3  (Phase 2)


 if t>=10
     va = 0;
 else
     va=20;
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% cell 1
  dcdt(1) = ((vsg1*ERK1^r/(Kag1^r+ERK1^r))+(vsg2*G1^s/(Kag2^s+G1^s)))*(Kig^q/(Kig^q+N1^q))-kdg*G1;

  dcdt(2) = ((vsn1*Kin1^u/(Kin1^u+ERK1^u))+(vsn2*N1^v/(Kan^v+N1^v)))*(Kin2^w/(Kin2^w+G1^w))-kdn*N1;

  dcdt(3) = vsfr1*Kifr/(Kifr+N1)+(vsfr2*G1)/(Kafr+G1)-kdfr*FR1;

  dcdt(4) = va*FR1*Fp1/(Kd+Fp1)*((1-ERK1)/(Ka+1-ERK1))-vi*ERK1/(Ki+ERK1);

  dcdt(5) = vex+(vsf*N1^z)/(Kaf^z+N1^z)-kdf*Fs1;

  

  end
