clear all;
clc;
%---------------------------------------------------------------%%
% PROJEKT TAP, KRZYSZTOF KASPRZYK, ADRIAN OPALACH, Artur Kropop %%
%                                                               %%
%---------------------------------------------------------------%%
% sta³e:
Ro = 10.^6 ;
Roc = Ro;
Cp = 1;
Cpc = 1;
Ko = 10.^10;
ER = 8330.1;
% ER = E/R
h = 130 * 10.^6;
a = 0.516*10.^6;
b = 0.5;
%--------------------------------------------
% Punkt pracy
V0 = 1 ;
Fin0 = 1;
F0 = 1;
Cain0 = 2;
Fc0 = 15;
Tin0 = 343;
Tcin0 = 310;
Ca0 = 0.16;
T0 = 405;
Ca_lin = Ca0;
T_lin = T0;
%--------------------------------------------------------------------------
dT = 0.1;
time = [0 : dT : 200]';

%Macierze wyjœciowe:

y_wyj = ones(size(time,1),2);
y_wyj(:,1) = Ca0*y_wyj(:,1);
y_wyj(:,2) = T0*y_wyj(:,2);

%wektor stanu pocz¹tkowego:

x0 = [Ca0 ; T0];

%wektor wejœciowy ( sterowanie):

u0 = [ Cain0 ; Fc0];

% skok sterowania:
u_step = [Cain0 ; Fc0];

%trajektoria zmian wejœæ obiektu:

u_traj= ones(size(time,1),2);
u_traj(:,1) = u_traj(:,1)*u0(1);
u_traj(:,2) = u_traj(:,2)*u0(2);

u_real_traj = ones(size(time,1),2);
u_real_traj(:,1) = u_traj(:,1)*u0(1);
u_real_traj(:,2) = u_traj(:,2)*u0(2);

for i=1:size(time,1)
    u_traj(i,1) = u_step(1);
    u_traj(i,2) = u_step(2);
    
    u_real_traj(i,1) = u_step(1);
    u_real_traj(i,2) = u_step(2);
    
end

%wektor zak³óceñ:

z0 = [Tin0 ; Tcin0];

% Skok zak³óceñ:

z0_step = [Tin0; Tcin0];

z_traj = ones(size(time,1),2);
z_traj(:,1)=z_traj(:,1)*z0(1);
z_traj(:,2)=z_traj(:,2)*z0(2);

for i = 1 : size(time,1)
    z_traj(i,1) = z0_step(1);
    z_traj(i,2) = z0_step(2);
end

u = u_traj(1,:)';
z = z_traj(1,:)';

%Macierze równan stanu:
global A;
global B;
global C;
A = zeros(2,2);
  A(1,1) = -12.6780;
  A(1,2) = -0.0949;
  A(2,1) = 1.5181e+03;
 A(2,2) = 7.4622; 
 % A(2,2) = 11.3359; 
  
B = zeros(2,4);
 B(1,1) = Fin0/V0;
 B(1,2) = 0;
 B(1,3) = 0;
 B(1,4) = 0;
 B(2,1) = 0;
 B(2,2) = -6.3038;
 B(2,3) = 1;
 B(2,4) = 1.8736; 
 %B(2,4) = 0.0664; 
C = zeros(2,2);
C(1,1)=1;
C(1,2)=0;
C(2,1) = 0;
C(2,2)=1;
D = zeros(2,1);

%Linearyzacja:
y_lin_traj = y_wyj;

x02 = x0;

x_lin = [0; 0];

u_lin = [ u(1)-u0(1);
    u(2)-u0(2);
    z(1)-z0(1);
    z(2)-z0(2);];
options = odeset('AbsTol',1e-6,'RelTol',1e-6);

y_lin = zeros(2,1);

time2= size(time,1);

sim_x0 = x0;
x_lin2 = zeros(2,1);
%g³ówna pêtla:
%func = [Fin0*u(1)/V0-F0*x(1)/V0-Ko*exp(-ER/x(2))*x(1) ; 
%    (1/(V0*Ro*Cp))*(Fin0*Ro*Cp*z(1)-F0*Ro*Cp*x(2)+V0*h*Ko*exp(-ER/x(2))*x(1)-a*u(2)^(b+1)/(u(2)+a*u(2)^b/(2*Roc*Cpc))*(x(2)-z(2)))];

for i=2 : time2
    
    sim_t=[time(i);time(i)+dT];
    [t , x]=ode45(@zbiornik,sim_t,sim_x0,options,u,z);
    sim_x0=x(size(x,1),:);
    u = u_traj(i,:)';
    z = z_traj(i,:)';
    
    y_wyj(i,1) = sim_x0(1);
    y_wyj(i,2) = sim_x0(2);
    
    [t , x_lin] = ode45(@linearyzacja,sim_t,x_lin,options,u_lin);
    
    x_lin = x_lin(size(x_lin,1),:);
    x_lin2(1)=x_lin(1);
    x_lin2(2) = x_lin(2);
    y_lin = C * x_lin2 ;
    y_lin_traj(i,:) = y_lin_traj(i,:) + y_lin'; 
    u_lin = [ u(1)-u0(1);
    u(2)-u0(2);
    z(1)-z0(1);
    z(2)-z0(2);];

end

%figure(1);
plot(time,y_wyj(:,2),'r','Linewidth',2);
hold on;
xlabel('Czas [s]');
plot(time,y_lin_traj(:,2),'b','Linewidth',2);
ylabel('Temperatura [K]');
legend('M. nieliniowy','M. liniowy');
axis([ 0 50 300 500]);
grid;
print('Wykres','-dpng','-r300');

figure(2);
plot(time,y_wyj(:,1),'r','Linewidth',2);
xlabel('Czas [s]' );
hold on;
plot(time,y_lin_traj(:,1),'b','Linewidth',2);
axis([0 50 -1 2]);
ylabel( 'Stê¿enie Ca [kmol/m3]');
legend('M. nieliniowy','M. liniowy');
grid;
print('Wykres2','-dpng','-r300');












