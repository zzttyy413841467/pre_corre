% clear;
constant_sim;

ee0=1/((Re+121900)/Re)-(7630/Vc)^2/2;
% eef=1/((Re+start(1))/Re)-(start(2)/Vc)^2/2;
eef=1/((Re+26500)/Re)-(900/Vc)^2/2;
espan=linspace(ee0,eef,5000);
yyy0=[55/57.3;(Re+121900)/Re;-0.5*pi/180;tau0;0;0;90/180*pi];
options=odeset('events',@stop_conditions);
[eee1,yyy1]=ode45(@dxde3,espan,yyy0,options);

e0=eee1(end);
ef=eef;
sigma0=sig0*pi/180;
sigmaf=40.5*pi/180;

espan=linspace(e0,ef,500);
s0=32/57.3;
% gamma0=start(3)*pi/180;
tau0=0;
y0=yyy1(end,:)';
% options=odeset('events',@stop_conditions);
[eee2,yyy2,sigma]=rk(@dxde1,espan,y0);

figure(1);
hold on
ee=[eee1;eee2];
yy=[yyy1;yyy2];

v=sqrt(2*(1./yy(:,2)-ee));
vvv=Vc*v;
hhh=Re*yy(:,2)-Re;
ttt=yy(:,4)*sqrt(Re/g0);
plot(vvv,hhh);
grid on

figure(2);
grid on
hold on
plot(Vc*v(1:end-1), [ones(length(eee1),1)*sig0;sigma*180/pi]);

zoulang;