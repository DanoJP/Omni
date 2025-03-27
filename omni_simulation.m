clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tf = 20; %Simulation Time
hs = 0.006; %Delta step Euler Method
h=hs;
o = Tf/hs;
t = 0:hs:Tf;
disp('Iniciando');
%% Parametros del Sistema %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mR=1.99;
mR1=0.29;
IRz=6.0848e-2;
IRy1=3.24e-4;
IRz1=4.69e-4;
r=0.05; 
Jm=5.7e-7;
kb=0.01336e-3; 
ka=0.0134;
Ra=1.9;
kv=0.0001; 
re= 64; 
L=.11;
E=-1/r*[sqrt(3)/2 -1/2 -L;0 1 -L;-sqrt(3)/2 -1/2 -L]';
B=-[0 -1 0;1 0 0;0 0 0];
mR11=mR+3*mR1;
mR33=3*mR1*L^2+IRz+3*IRz1;
MR=[mR11 0 0;0 mR11 0;0 0 mR33];
M=MR+(IRy1+Jm*re^2)*(E'*E);
D=re^2*(ka*kb/Ra+kv)*(E'*E);
G=inv(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Trayectoria deseada
w=2*pi*(0.05);
w2=0.5*pi*(0.05);

xd=0.5*sin(w*t);
yd=0.5*cos(w*t);
thetad= -w2*t;
    
xpd=0.5*w*cos(w*t);
ypd=-0.5*w*sin(w*t);
thetapd=-w2+0*t;% se agrega 0*t para que sea consistente respecto al tiempo

xppd=-0.5*w*w*sin(w*t);
yppd=-0.5*w*w*cos(w*t);
thetappd=0*t;

xid = [xd ; yd ; thetad] ;
xipd = [xpd ; ypd ; thetapd] ;
xippd = [xppd ; yppd ; thetappd] ;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Disturbances
delta(1,:) = 0.01*sin((pi/10)*t) +0.01;
delta(2,:) = 0.01*cos((pi/10)*t) +0.01;
delta(3,:) = 0.01*cos((pi/10)*t) +0.01;
%% Bound of Disturbances
eta1 = .02*(pi/10)*5;
eta2 = .02*(pi/10)*5; 
eta3 = .02*(pi/10)*5; 

%% Controler Gains
    ceta= 1.1;
    etab(1,1) = ceta*eta1;
    etab(2,1) = ceta*eta2;
    etab(3,1) = ceta*eta3;
    %% NTSMC Gains
        c(:,1)   = 5;
        kb1(:,1) = 16.*etab.^(2/3)*1;
        kb2(:,1) =   7.*etab*1;
    
    %% Observer gains
    k(1,1) = 2*eta1^(1/3);
    k(2,1) = 2*eta2^(1/3);
    k(3,1) = 2*eta3^(1/3);

    k3(1,1) = 1.5*eta1^(1/2);
    k3(2,1) = 1.5*eta2^(1/2);
    k3(3,1) = 1.5*eta3^(1/2);

    kz(1,1) = 1.1*eta1;
    kz(2,1) = 1.1*eta2;
    kz(3,1) = 1.1*eta3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Condiciones Iniciales
x(1)= 0;
x_p(1)=0;
x_pp(1)=0;

y(1)= 0;
y_p(1)=0;
y_pp(1)=0;

theta(1)=0;
theta_p(1)=0;
theta_pp(1)=0;

xi1(:,1)=[x(1);y(1);theta(1)];
xi2(:,1)=[x_p(1);y_p(1);theta_p(1)];

 %R(:,:,1) = [0 0 0; 0 0 0; 0 0 0];

xh(:,1) = [0;0;0];
xph(:,1) =  [0;0;0];
% xp(:,1) =  [0;0;0];
sh(:,1) = [0;0;0];
zeta(:,1) = [0;0;0];

v(:,1) = [0;0;0];


    for i=1:o
    %% Tracking errors
    e(:,i) = xid(:,i)- xi1(:,i); % epsilon pi
    eph(:,i) = xipd(:,i)-xph(:,i);% epsilon vi
    ee(:,i) = xi1(:,i)- xh(:,i); %epi
    % epe(:,i) = xp(:,i)-xph(:,i);
    %% Matrix calculations
    R(:,:,i)=[cos(xi1(3,i)) sin(xi1(3,i)) 0;-sin(xi1(3,i)) cos(xi1(3,i)) 0;0 0 1];
    C(:,:,i)=2/(r^2)*(IRy1+Jm*re^2)*xph(3,i)*B;
    F(:,i)=-inv(M)*((C(:,:,i)+D)*[xph(1,i),xph(2,i),xph(3,i)]');
    
    c_d     = 0;  % Constant pertubation 


sh(:,i) =e(:,i) + c.*(abs(eph(:,i)).^(3/2).*sign(eph(:,i)));
v(:,i+1) = v(:,i) + (-kb2.*(abs(sh(:,i)).^(0).*sign(sh(:,i))))*h;
u2(:,i) =k3.*(abs(ee(:,i)).^(1/3).*sign(ee(:,i))) - kb1.*(abs(sh(:,i)).^(1/3).*sign(sh(:,i))) + v(:,i); 
    
tau(1,i)= -(1/G(1,1))*(F(1,i)-xippd(1,i) + zeta(1,i) + u2(1,i) + c_d);
tau(2,i)= -(1/G(2,2))*(F(2,i)-xippd(2,i) + zeta(2,i) + u2(2,i) + c_d);
tau(3,i)= -(1/G(3,3))*(F(3,i)-xippd(3,i) + zeta(3,i) + u2(3,i) + c_d);

%% Observer 
    xh(:,i+1) = xh(:,i) + (xph(:,i) + k.*abs(ee(:,i)).^(2/3).*sign(ee(:,i)))*h;
    xph(:,i+1) = xph(:,i) + (zeta(:,i) + k3.*(abs(ee(:,i)).^(1/3)).*sign(ee(:,i)) +F(:,i) + diag(G).*tau(:,i))*h;
    zeta(:,i+1) = zeta(:,i) + (kz.*((abs(ee(:,i)).^(0)).*sign(ee(:,i))))*h;
    %% Dynamical model
    xi1(:,i+1)= xi1(:,i)+xi2(:,i)*h;
    % xp(1,i+1)= xp(1,i) +(F(1,i)+G(1,1)*(delta(1,i) + tau(1,i)))*h;
    % xp(2,i+1)= xp(2,i) +(F(2,i)+G(2,2)*(delta(2,i) + tau(2,i)))*h;
    % xp(3,i+1)= xp(3,i) +(F(3,i)+G(3,3)*(delta(3,i) + tau(3,i)))*h;
    xi2(:,i+1) = (F(:,i) +diag(G).*(delta(:,i) + tau(:,i))) * hs +xi2(:,i);
modulo = mod(i,100);    
if modulo==0  

        %% Graph
figure(2)
clf
plot(xid(1,1:i),xid(2,1:i),'k','linewidth',1)
hold on
plot(xi1(1,1:i),xi1(2,1:i),'-.','linewidth',1.5)
hold on
% legend('x_d','CNTSMC','Fontsize',10)
xlabel('x_1')
ylabel('x_2')
axis([-2.5 2.5 -2.5 2.5]);
% axis([-1 1 -1 1]);
end

    end

%% Grafica X_1
subplot(1,3,1)
plot(t,xid(1,:),'k','linewidth',1)
hold on
plot(t,xi1(1,:),'-.','linewidth',1.5)
hold on
ylabel('Metros')
xlabel('Tiempo[s]')
title('x')

%% Grafica X_2
subplot(1,3,2)
plot(t,xid(2,:),'k','linewidth',1)
hold on
plot(t,xi1(2,:),'-.','linewidth',1.5)
hold on
ylabel('Metros')
xlabel('Tiempo[s]')
title('y')


%% Grafica X_3
subplot(1,3,3)
plot(t,rad2deg(xid(3,:)),'k','linewidth',1)
hold on
plot(t,rad2deg(xi1(3,:)),'-.','linewidth',1.5)
hold on
title('\theta')
xlabel('Tiempo[s]')
ylabel('Grados')
legend({'x_d','CNTSMC'},'Orientation','horizontal','Fontsize',10)
