function [value,isterminal,direction]=stop_conditions(e,x)

constant_sim;
r=x(1);
gamma=x(4);
v=sqrt(2*(1./r-e));
alpha0=45/180*pi;
ma=v*Vc/340;
alpha=alpha0.*(ma>=15)+((45-0.21*(ma-15).^2)/180*pi).*(ma<15);
Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;
h=r*Re-Re;
rho=rho0.*exp(-h./hs);
q=1/2*rho.*v.^2;
L=q.*Cl.*S/m*Re;
D=q.*Cd.*S/m*Re;

drdv=v.*sin(gamma)./(-D-sin(gamma)./r.^2);

f=@(r1)1./r1.^2-v.^2./r1-1/2*rho0.*exp(-(r1*Re-Re)./hs).*v.^2.*Cl.*S./m.*Re;
r1=fsolve(f,1,optimset('Display','off'));
K=Re*S*Cl/2/m;
rho1=rho0.*exp(-(r1*Re-Re)./hs);
drhodr=-Re/hs*rho1;
% drdvq=-2*v.*(K.*rho1.*r1.^3+r1.^2)./(K.*r1.^3.*drhodr.*v.^2-v.^2.*r1+2);
%drdvq=-2*v.*(K.*rho1.*r1.^2+r1).^2./(K.*r1.^2.*drhodr+2*K.*rho1.*r1+1);
% drdvq=-6.3*q_max^2/(k_q^2*drhodr*v^7.3);
drdvq=-6.3*rho1/(drhodr*v);
delta=0.01;

value=abs(drdv-drdvq)-delta;
isterminal=1;
direction=-1;
end
