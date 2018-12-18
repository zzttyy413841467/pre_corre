clear;
constant_sim;
calcu_sigma0;

thetaf=65*pi/180;
phif=5*pi/180;

global nfz;
nfz=struct('theta',{25/180*pi 45/180*pi},'phi',{4/180*pi 0.8/180*pi},'radius',{2/180*pi 2.5/180*pi});

ee0=1/((Re+121900)/Re)-(7630/Vc)^2/2;
eef=1/((Re+26500)/Re)-(900/Vc)^2/2;

r0=(Re+121900)/Re;
theta0=0;
phi0=0;
gamma0=-0.5*pi/180;
psi0=88/180*pi;
tau0=0;
% psi_t=atan(sin(thetaf-theta0)/(cos(phi0)*tan(phif)-sin(phi0)*cos(thetaf-theta0)));
if (phif-phi0>0)
    psi_t=asin(cos(phif)*sin(thetaf-theta0)/sin(acos(cos(phi0)*cos(phif)*cos(theta0-thetaf)+sin(phi0)*sin(phif))));
else
    psi_t=pi-asin(cos(phif)*sin(thetaf-theta0)/sin(acos(cos(phi0)*cos(phif)*cos(theta0-thetaf)+sin(phi0)*sin(phif))));
end

sig0=-sig0_a*sign(psi0-psi_t);
sig_0=sig0;

yyy0=[r0;theta0;phi0;gamma0;psi0;tau0];
options=odeset('events',@stop_conditions);
[eee1,yyy1]=ode45(@dxde3,linspace(ee0,eef,2000),yyy0,options);
% zoulang;
% figure(1)
% v1=sqrt(2*(1./yyy1(:,1)-eee1));
% vvv1=Vc*v1;
% hhh1=Re*yyy1(:,1)-Re;
% plot(vvv1,hhh1);

e0=eee1(end);
ef=eef;
sigma0=sig0_a*pi/180;
sigmaf=40*pi/180;
espan=linspace(e0,ef,300);
y0=yyy1(end,:)';
[eee2,yyy2,sigma]=rk(@dxde1,espan,y0);

figure(1);
hold on
ee=[eee1;eee2];
yy=[yyy1;yyy2];

vv=sqrt(2*(1./yy(:,1)-ee));
vvv=Vc*vv;
hhh=Re*yy(:,1)-Re;
ttt=yy(:,6)*sqrt(Re/g0);
plot(vvv,hhh);
grid on

figure(2);
grid on
hold on
plot(Vc*vv(1:end-1), [ones(length(eee1),1)*sig_0;sigma*180/pi]);

zoulang;

figure(3);
hold on
theta=0:1/50*pi:2*pi;
x1=nfz(1).theta+nfz(1).radius.*cos(theta);
y1=nfz(1).phi+nfz(1).radius.*sin(theta);
x2=nfz(2).theta+nfz(2).radius.*cos(theta);
y2=nfz(2).phi+nfz(2).radius.*sin(theta);
plot(x1*180/pi,y1*180/pi);
plot(x2*180/pi,y2*180/pi);
plot(180/pi*yy(:,2),180/pi*yy(:,3))
grid on

% openfig('nfz.fig');
