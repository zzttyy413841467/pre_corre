function [e,x1,sigma]=rk(dxde1,espan,x0)
n=length(espan);

ef=espan(end);
e0=espan(1);
h=(ef-e0)/(n-1);
i=1;
x=x0;
global e_m
global sigmam
e_m=(ef-e0)/5+e0;

x1=zeros(n,length(x0));
sigmaf=evalin('base','sigmaf');
sigma=zeros(n-1,1);
x1(1,:)=x;
sigma0=evalin('base','sigma0');
sigmam=sigma0;
sf=zeros(n,1);
sigma_m=zeros(n-1,1);
sigma_m(1)=sigmam;
for e=espan(1:end-1)
    if i==2
        sigmam=sigmam-2/57.3;
        sigma_m(i)=sigmam;
    end
    if i>2
        sigmam=magnitude(sigma_m(i-1)-(sigma_m(i-1)-sigma_m(i-2))/(sf(i-1)-sf(i-2))*sf(i-1));
        sigma_m(i)=sigmam;
    end
    [ee1,xx1]=ode45(@dyde,linspace(e,ef,200),x(1:4));
    sf(i)=xx1(end,1);
    if(e<=e_m)
        sigma_x=limit(sigma0+(sigmam-sigma0)./(e_m-e0).*(e-e0),e,x(1:4));
    else
        sigma_x=limit(sigmam+(sigmaf-sigmam)./(ef-e_m).*(e-e_m),e,x(1:4));
    end
    sigma(i)=sigma_x;
    k1=dxde1(e,x);
    k2=dxde1(e+h/2,x+h*k1/2);
    k3=dxde1(e+h/2,x+h*k2/2);
    k4=dxde1(e+h,x+h*k3);
    x=x+(h/6)*(k1+2*k2+2*k3+k4);
    i=i+1;
    x1(i,:)=x;

end
e=espan';
% x1(((i+1):n),:)=[];
end