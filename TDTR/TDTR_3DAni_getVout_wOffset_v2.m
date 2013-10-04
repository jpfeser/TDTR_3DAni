function Vout = TDTR_3DAni_getVout_wOffset_v2(tdelay,k,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%require tdelay to be a column vector
tdelay = tdelay(:);

fmax=10/min(abs(tdelay)); %maximum frequency considered (see RSI paper)
ii=sqrt(-1);

M=10*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
mvect=-M:M; %Range of Fourier components to consider (Vector)
fudge1=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
fudge2=exp(-pi*((mvect/tau_rep-f)/fmax).^2);

dT1=zeros(1,length(mvect))';
dT2=zeros(1,length(mvect))';
boundary_x = 1/sqrt(ws_x^2+wp_x^2)*2;
boundary_y = 1/sqrt(ws_y^2+wp_y^2)*2;
a=-2*pi*boundary_x;
b=2*pi*boundary_x;
c=-2*pi*boundary_y;
d=2*pi*boundary_y;

fprintf('calculating freq responses...\n')
tic
if license('test','distrib_computing_toolbox')
parfor i = 1:length(mvect)
    [dT1(i),n1]=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral_wOffset_v2(ETA,NU,FF,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset),a,b,c,d,mvect(i)/tau_rep+f);
    [dT2(i),n2]=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral_wOffset_v2(ETA,NU,FF,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset),a,b,c,d,mvect(i)/tau_rep-f);
    %[i,2^(2*n1)] %comment to hide progress output
end
else
for i = 1:length(mvect)
    [dT1(i),n1]=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral_wOffset_v2(ETA,NU,FF,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset),a,b,c,d,mvect(i)/tau_rep+f);
    [dT2(i),n2]=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral_wOffset_v2(ETA,NU,FF,k,C,h,wp_x,wp_y,Qp,ws_x,ws_y,xoffset,yoffset),a,b,c,d,mvect(i)/tau_rep-f);
    %[i,2^(2*n1)] %comment to hide progress output
end
end
toc
fprintf('finished freq responses...\n')

dT1 = conj(dT1); % I don't know why I should have to do this?
dT2 = conj(dT2);

expterm=exp(ii*2*pi/tau_rep*(tdelay*mvect));
Retemp=(ones(length(tdelay),1)*(dT1'.*fudge1+dT2'.*fudge2)).*expterm;
Imtemp=-ii*(ones(length(tdelay),1)*(dT1-dT2)').*expterm;

Resum=sum(Retemp,2); %Sum over all Fourier series components
Imsum=sum(Imtemp,2);

deltaRm=TCR*(Resum+ii*Imsum); %
deltaR=deltaRm.*exp(ii*2*pi*f*tdelay); %Reflectance Fluxation (Complex)

ratio=-real(deltaR)./imag(deltaR);
Vout = imag(deltaR);

end

