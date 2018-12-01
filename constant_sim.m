m=35828;
S=149.4;

q_max=790000;

cl1 = 0.016292*(180/pi);                  % Parameters for lift coefficient
cl2 = 0.00026024*(180/pi)^2;       
cl0 = -0.041065;
cd1 = -0.03026;       
cd2 = 0.86495;       
cd0 = 0.080505;

rho0=1.225;
g0=9.81;
Re=6371393;
mu=g0*Re^2;
hs=7320;
Vc=sqrt(g0*Re);

k_q=9.4369e-5*Vc^3.15;