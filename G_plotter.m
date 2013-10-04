%%Generate materials properties
clear all

%% More Complicated
lambda = [142 1420];
phi = 1*ones(size(lambda));
lambda_tensor = MakeRadial_ktensor(lambda,phi);
C = [1.62 1.62]*1e6;
h = [100 1e6]*1e-9;
%lambda_tensor(:,1:2)=0;
TCR = 1e-4;

%% Generate laser parameters
wp = 2.5e-6;
ws = 2.5e-6;
Qp = 30e-3*pi/2;

%% Generate integration parameters
kmax=1/sqrt(wp^2+ws^2)*2;
a = -2*pi*kmax;
b = 2*pi*kmax;
c = -2*pi*kmax;
d = 2*pi*kmax;

%%Generate the spacings to plot
f = 10e3;
eta = 0;
xi = logspace(log10(b/1e4),log10(b),200);
kvect = xi/(2*pi);
[XI,ETA,FF]=ndgrid(xi,eta,f);

%% Evaluate G, P, S
[G] = TDTR_3DAni_getG_savespace(XI,ETA,FF,lambda_tensor,C,h);
P = 1/(2*pi)*exp(-1/8*wp^2*(xi.^2+eta.^2)); 
S = 1/(2*pi)*exp(-1/8*ws^2*((-xi).^2+(-eta).^2));

%% Ditto for non-3D code
[Integrand,G_radial]=TDTR_TEMP_VV3(kvect,f,lambda,C,h,phi,wp,ws,Qp);

%% plot shape of result
hold off
semilogx(XI(:),G(:)/max(G(:)),'g-');
hold on
semilogx(XI(:),P(:)/max(P(:)),'b-');
semilogx(XI(:),S(:)/max(S(:)),'r-');
semilogx(2*pi*kvect(:),G_radial(:)/max(G_radial(:)),'go');
axis square
hold off

%% Compare 3D, Radial Solutions
hold off
semilogx(XI(:),G(:),'g-');
hold on
semilogx(2*pi*kvect(:),G_radial(:),'go');
axis square
hold off