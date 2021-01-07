%% wind turbine blade design solver
%% Authors
%  Yuvarajendra Anjaneya Reddy(yuvan983)
%  Rohith Prem Maben(rohma417)


%%                                 Aerodynamics assignment
%% First submission  15 /10 /2018

%% Updated for resubmission 10/31/2018
%   (1)changes chord and pitch distribution from ideal case
%   (2)added ideal and actual data for graphs for comparison

tic
clear all
close all
clc

Rtip=45;                         %Radius of Rotor
B=3;                             %Number of blades
U=10                            %wind velocity
rho=1.225;                       %Fluid density
N=100;                           %Number of blade elements
tolerance=0.0001;                %Tolerance for iteration
dr=Rtip/N;                       %Element width
Vrange=5:1:26;                   %off design velocity condition

Rhub=0.05*Rtip;                %5 percent ignored for hub location

r=linspace(Rhub,Rtip,100);
Rr=r/Rtip         ;            %non dimensional blade radius

               

%% Blade element pitch and chord distribution 

%% ideal shape assumption
omega_i=1.884;                       % rad/sec
lam_est_i=(Rtip*omega_i)/U;          


Cl_root_design=1.4;
aoa_root_design=8;

Cl_tip_design=0.841;
aoa_tip_design=6;

for i=1:100
lam_i(i)=lam_est_i*(r(i)/Rtip);

phi_design(i)=(2/3).*atan(1/lam_i(i));
phi_deg_design(i)=rad2deg(phi_design(i));
 if i<50

   C_1(i)=(8*pi*r(i)*(1-cos(phi_design(i))))./(B.*Cl_root_design);
   aoa_ideal(i)=aoa_root_design;
   Re_design(i)=2.5*10.^6;
 elseif i>=50
   C_1(i)=(8*pi*r(i)*(1-cos(phi_design(i))))./(B.*Cl_tip_design);
   aoa_ideal(i)=aoa_tip_design;
   Re_design(i)=3*10.^6;
 end
end

% linear fit to ideal chord distribution
p = polyfit(r,C_1,1); 
a1=p(1);b1=p(2);
chord_est=a1*r+b1;
hold on
% quadratic fit to ideal chord distribution
p2 = polyfit(r,C_1,2); 
a2=p2(1);b2=p2(2);
c=p2(3);
chord_est_2=a2*r.^2+b2.*r+c;

theta0=0;
thetai=phi_deg_design-aoa_root_design;
theta_ideal=thetai-theta0;

%%
pitch_dist     = theta_ideal;
pitch_dist_rad = pitch_dist.*(pi/180);


C=chord_est;

%% store 2 data files for airfoils

% S818 design Re 2.5 mil
root_aoa=xlsread('airfoildata.xlsx','A4:A22');
root_Cd =xlsread('airfoildata.xlsx','C4:C22');
root_Cl =xlsread('airfoildata.xlsx','B4:B22');
% S828 design Re 3 mil
tip_aoa=xlsread('airfoildata.xlsx','E4:E15');
tip_Cd =xlsread('airfoildata.xlsx','G4:G15');
tip_Cl =xlsread('airfoildata.xlsx','F4:F15');

Drag_r=root_Cd;
Lift_r=root_Cl;
root_ld=Lift_r.*Drag_r.^-1;

Drag_t=tip_Cd;
Lift_t=tip_Cl;
tip_ld=Lift_t.*Drag_t.^-1;
%% angular velocity and tip speed ratio

%omega=1.884              % rad/sec
lambda=7;                 % tip speed assumed
lambda_Ratio=lambda*(Rr);
omega =2.6;              %rad/sec angular velocity
%% initial assumption for a and a'
w           = U*lambda*Rtip.^-1;
a(1,1:100)   = 0.333;%0.333;
angular(1,1:100)=(1-(3*a(1))).*((4*a(1))-1).^-1;
 for i=1:length(r)
      lambda_section(i)=(w*r(i))*U.^-1;
      phi_int(i)       =(2*3.^-1)*(atan(1*(lambda_section(i).^-1)).*180*pi.^-1);
      
 end

 
 % error initialize
delta_axial  =1;
delta_angular=1;

itr=0;


%% iteration for 15 elements along the span

for N=1:100
    
   while delta_axial> tolerance & delta_angular > tolerance
    
       itr=itr+1;
      
       old_a   = a;
       old_aa  = angular;
%  1:  angle of relative wind
phi              = atan((1-a)./((1+angular ).*(lambda_Ratio ))); % eq 3.5.6
phi_degrees      = phi .*(180/pi); 


%  2:  Tip loss correction
F =(2/pi)*acos(exp(-(((B/2).*(Rtip-(r)))./((r).*sin(phi)))));    %eq 3.7.19


%  3:  angle of attack
alpha=phi-pitch_dist_rad;                          %eq3.5.5
alpha_degrees=alpha.*(180/pi);


%% Calculate Lift and Drag coefficients for angle of attack >20 degrees
            if N<=7                               % root section

                    if  alpha_degrees >20
                           C_d=2.*sin(phi_degrees).^2;                 %flat plate airfoil
                           C_l=2.*sin(phi_degrees).*cos(phi_degrees);  % flat plate airfoil
                           
                      else C_d=interp1(root_aoa,Drag_r,(alpha));
                           C_l=interp1(root_aoa,Lift_r,(alpha));
                    end
           
             elseif N>7 %tip section

                    if alpha_degrees >20
                           C_d=fprintf('aoa outside blade envelope\n');%2.*sin(phi_degrees).^2;
                           C_l=fprintf('aoa outside blade envelope\n');%2.*sin(phi_degrees).*cos(phi);
                           
                      else C_d=interp1(tip_aoa,Drag_t,(alpha));
                           C_l=interp1(tip_aoa,Lift_t,(alpha));
                    end
            end
           
           
% force and torque coefficients
C_n=(C_l.*cos(phi))+(C_d.*(sin(phi)));
C_t=(C_l.*sin(phi))-(C_d.*(cos(phi)));
 

solidity=(B.*C)/(2*pi.*Rtip);

% 6: Calculate new values for induction factors
a       =((1+((4.*F(N)*sin(phi).^2).*(solidity.*C_n).^-1))).^-1;
angular =1./(((4.*F(N).*sin(phi).*cos(phi))./( solidity.*C_t))-1);

% 7:update old values
a_new  = a;
aa_new = angular;

% 8: Calculate the difference between iterations
delta_axial=abs(old_a-a_new);
delta_angular=abs(old_aa-aa_new);
      
    end   



end
dT=4*pi.*r.*rho*U.^2.*a_new.*(1-a_new).*dr.*F;       %Thrust (3.5.1)
dQ=4*pi.*r.^3.*rho.*U.*omega.*aa_new.*(1-a_new).*dr.*F; %Torque (3.5.2)

dP=dQ*omega;                                         %Power
M=sum(dQ);
P=sum(dP)*10.^-6;                                    % power output in megawatts
P_a=0.5*rho*pi*Rtip^2*U^3*10.^-6;                    % power available in megawatts
C_p=P/P_a;                                           % power coefficient
Urelative=U.*(1-a_new).*(sin(phi)).^-1;
%% plotting of results

%% plot 1 chord distribution
figure
plot(r,chord_est,r,chord_est_2,r,C_1)
title('chord distribution along radius')
xlabel('radius(m)')
ylabel('chord distribution(m)')
legend('linear fit','ideal case','quadratic fit') 
print('chord', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
%% plot 3 blade veloicties
figure
% f_high=(12)/60  ;% 1/s% from literature of similar category turbine(vestasV90)
% w_rps=2*pi*f_high;
 v_angular=omega*r;%*(1/(2*pi));%r*w_rps;
Urel=sqrt(U.^2+(v_angular).^2);
Udesign=10*ones(1,100);
plot(r,v_angular,r,Urelative,r,Udesign)
title('angular veloicty vs blade radius')
xlabel('radius(m)')
ylabel('velocity(m/s)')
legend('angular velocity(m/s)','relative wind veloicty','Design wind veloicty')
print('velocities', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

%% plot 2 reynolds number distribution
figure 
mu=1.785*10.^-5;
Re_linear=(rho*U.*chord_est)*(mu.^-1);
%(rho*U*chord_est)*(mu.^-1);

plot(r,Re_design);
ylim([1*10.^6  4*10.^6])
title('Reynolds number distribution along radius')
xlabel('radius(m)')
ylabel('reynolds number')
legend('Reynolds number(simplified)')


print('reynolds', '-dpng', '-r300'); %<-Save as PNG with 300 DPI



%% 2a blade aoa vs radius
figure 
plot(r,alpha_degrees,r,pitch_dist,r,phi_degrees,r,aoa_ideal)
title('Blade angles')
xlabel('radius(m)')
ylabel('sectional angle of attack(degree)')
legend('angle of attack','pitch angle','relative angle','ideal angle of attack')
print('blade angles', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

figure 
plot(r,aoa_ideal,r,alpha_degrees)
title('Blade angle of attack(deg)')
xlabel('radius(m)')
ylabel('sectional blade pitch angle(degree)')
legend('calculated angle of attack','design angle of attack')


%% 2b blade pitch vs radius
figure 
plot(r,pitch_dist)
title('Blade pitch distribution(deg)')
xlabel('radius(m)')
ylabel('sectional blade pitch angle(degree)')
print('pitch', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
%% 2c relative wind angle vs radius
figure 
plot(r,phi_degrees,r,phi_deg_design)
title('Blade phi distribution(deg)')
xlabel('radius(m)')
ylabel('relative wind angle')
legend('relative wind angle(calculated)','relative wind angle(ideal)')
print('relative wind angle', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
%% 2d  L D ratio vs aoa
 figure
 title('lift over drag ratio vs angle of attack')
 LDmax_root=max(root_ld);
 LDmax_tip=max(tip_ld);
 subplot(1,2,1)
 plot(root_aoa,root_ld)
 xlabel('angle of attack') 
 ylabel('lift to drag ratio')
 title('root airfoil:S818')
 subplot(1,2,2)
 plot(tip_aoa,tip_ld)
  xlabel('angle of attack') 
 ylabel('lift to drag ratio')
 title('tip airfoil:S828')
 
 print('aoa LD', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
% 
 %% 4a lift and drag coeff
figure 
plot(r,C_l,r,C_d)
legend('CL','CD')
xlabel('radius(m)')
ylabel('lift and drag coefficient')
title(' lift and drag coefficient')
 print('lift and drag coeff', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
% 
% 
% 4b force coeff
figure 
title('force coefficient')
plot(r,C_t,r,C_n)
legend('C tangential','C normal')

xlabel('radius(m)')
ylabel('tangential/normal force coefficient')
title('force coefficient')
print('force coefficients', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
% 4c thrust and force coeff
figure
plot(r,dT,r,dQ)
legend('torque','thrust')
title('blade forces')
xlabel('radius(m)')
ylabel('torque/thrust vs. blade radius')
print('forces', '-dpng', '-r300'); %<-Save as PNG with 300 DPI
% axial and angular induction factor
figure 
plot(r,a,r,angular)
xlabel('radius(m)')
ylabel('axial,angular induction factor')
title('axial,angular induction factor')
print('induction factors', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

% 5 tip loss factor
figure
plot(r,F)

title('tip loss correction factor')
xlabel('radius(m)')
ylabel('tip loss correction factor')
print('tiploss', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

%% 6 cp variation
figure
% U is varied to off design speeds and the following vector is created
Cp_vrange=[0.510 0.509 0.504 0.4978 0.485 0.4774 0.4646 0.4506 0.4488 0.4358 0.4267 0.4028 0.3779 0.3523 0.3267 0.3015 0.2772 0.2541 0.2324 0.2122 0.1936 0.1766]
u=5:1:26;  % off design speeds
%P_avail=0.5*rho.*pi.*Rtip.^2*u.^3*10.^-6;% power available in megawatts
%C_p_1=P.*P_avail.^-1;
plot(u,Cp_vrange,'-*')

title('variation of power coefficient with wind speed(off design)')
xlabel('windspeed(m/s)')
ylabel('power coefficient')
print('powercoeff', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

toc
