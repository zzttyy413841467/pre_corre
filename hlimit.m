function hmin=hlimit(e,x)

constant_sim;

r=x(2);
v=sqrt(2*(1/r-e));

alpha0=45/180*pi;
ma=v*Vc/340;
alpha=alpha0.*(ma>=15)+((45-0.21*(ma-15).^2)/180*pi).*(ma<15);

Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;

Kl=Re*Cl*S/2/m;

hqdot=2*hs.*log(sqrt(1.225).*(Vc*v).^3.15*9.4369e-5/750000);
hq=hs*log(1.225*(Vc*v).^2/2/1.5e4);
hn1=hs*log(1.225*S*(Vc*v).^2.*(Cl.*cos(alpha)+Cd.*sin(alpha))./2/2.5/m/g0);

hmin=max([hqdot hq hn1]);


end