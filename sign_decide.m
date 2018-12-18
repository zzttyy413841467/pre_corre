function y = sign_decide(sign0,e,x)

global nfz;
constant_sim;

delt_sigma1=10/180*pi;
delt_sigma2=20/180*pi;
delt_sigma3=5/180*pi;

r=x(1);
theta=x(2);
phi=x(3);
v=sqrt(2*(1./r-e))*Vc;
psi=x(5);
thetaf=evalin('base','thetaf');
phif=evalin('base','phif');

if(phif-phi>0)
    psi_t=asin(cos(phif)*sin(thetaf-theta)/sin(acos(cos(phi)*cos(phif)*cos(theta-thetaf)+sin(phi)*sin(phif))));
else
    psi_t=pi-asin(cos(phif)*sin(thetaf-theta)/sin(acos(cos(phi)*cos(phif)*cos(theta-thetaf)+sin(phi)*sin(phif))));
end
delt_sigma=psi-psi_t;

flag=0;
n=length(nfz);
for i=1:n
   if (nfz(i).theta>theta && sqrt((nfz(i).theta-theta)^2+(nfz(i).phi-phi)^2)-2*nfz(i).radius<0)
       flag=i;
       break;
   end  
end

if (flag==0)
    if (v>=5000)
        if (delt_sigma>delt_sigma1)
            sign=-1;
        elseif (delt_sigma<=delt_sigma1)&&(delt_sigma>=-delt_sigma1)
            sign=sign0;
        else 
            sign=1;
        end
    elseif (4000<=v&&v<5000)
        if (delt_sigma>delt_sigma2)
            sign=-1;
        elseif (delt_sigma<=delt_sigma2)&&(delt_sigma>=-delt_sigma2)
            sign=sign0;
        else 
            sign=1;
        end
    elseif (3000<=v&&v<4000)
        delt=delt_sigma3+(delt_sigma2-delt_sigma3)/(4000-3000)*(v-3000);
        if (delt_sigma>delt)
            sign=-1;
        elseif (delt_sigma<=delt)&&(delt_sigma>=-delt)
            sign=sign0;
        else 
            sign=1;
        end  
    else
        if (delt_sigma>delt_sigma3)
            sign=-1;
        elseif (delt_sigma<=delt_sigma3)&&(delt_sigma>=-delt_sigma3)
            sign=sign0;
        else 
            sign=1;
        end  
    end
else
    psi_1=asin(cos(nfz(flag).phi)*sin(nfz(flag).theta-theta)/sin(acos(cos(phi)*cos(nfz(flag).phi)*cos(theta-nfz(flag).theta)+sin(phi)*sin(nfz(flag).phi))));
    gamma=asin(nfz(flag).radius/sqrt((nfz(flag).theta-theta)^2+(nfz(flag).phi-phi)^2));
    if(psi_t>psi_1)
        if (psi<psi_1+gamma)
            sign=1;
        elseif (psi<=psi_1+30/180*pi+gamma && psi>=psi_1+gamma)
            sign=sign0;
        else
            sign=-1;
        end
    else
        if (psi<psi_1-30/180*pi-gamma)
            sign=1;
        elseif (psi>=psi_1-30/180*pi-gamma && psi<=psi_1-gamma)
            sign=sign0;
        else
            sign=-1;
        end
    end
    
end

if(sign0==0)
    sign=1;
end
y =sign;
end