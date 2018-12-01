constant_sim;
V=(850:50:8000)';
alpha0=45/180*pi;
ma=V/340;
alpha=alpha0.*(ma>=10)+((45-0.612*(ma-10).^2)/180*pi).*(ma<10);
Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;

hqdot=2*hs.*log(sqrt(1.225).*V.^3.15*9.4369e-5/q_max);
hq=hs*log(1.225*V.^2/2/1.5e4);
hn=hs*log(1.225*S*V.^2.*sqrt(Cl.^2+Cd.^2)./2/2.5/m/g0);
hn1=hs*log(1.225*S*V.^2.*(Cl.*cos(alpha)+Cd.*sin(alpha))./2/2.5/m/g0);
f=@(x)1/2*1.225*exp(-(x-Re)/hs).*V.^2.*Cl.*S./m*cos(10*pi/180)-mu./x.^2+(V.^2./x);
hb=fsolve(f,Re*ones(length(V),1),optimset('Display','off'))-Re;

figure(1);
plot(V,hqdot,':');
hold on
plot(V,hq,'-.');
hold on
plot(V,hn1,'--');
hold on
plot(V,hb);
grid on
ylim([0 1.2e5]);
ylabel('高度/m');
xlabel('速度/(m/s)');
legend({'热流边界','动压边界','过载边界','平衡滑翔边界'},'Location','NorthWest','FontSize',10);



v=V/Vc;
Kl=Re*Cl*S/2/m;

sigmaqdot=180/pi*acos(k_q.^2.*(1-v.^2).*v.^4.3./(Kl*q_max^2));
sigmaq=180/pi*acos(Vc^2*(1-v.^2)./(2*Kl*1.5e4));
sigman1=180/pi*acos((1-v.^2)./2.5.*(cos(alpha)+Cd./Cl.*sin(alpha)));

sigma_max1=min(min(sigmaqdot,sigmaq),sigman1);
figure(2);

plot(V,sigma_max1,'--');
hold on
plot(V,10*ones(1,length(V)));
grid on
ylim([0 90]);
ylabel('sigma/degree');
xlabel('速度/(m/s)');
legend({'热流边界','动压边界','过载边界','平衡滑翔边界'},'Location','NorthWest','FontSize',10);

