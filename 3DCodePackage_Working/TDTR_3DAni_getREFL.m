function [ratio,deltaR] = TDTR_3DAni_getREFL(tdelay,k,C,h,f,tau_rep,wp,Qp,ws,TCR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fmax=10/min(abs(tdelay)); %maximum frequency considered (see RSI paper)
ii=sqrt(-1);

M=10*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
mvect=-M:M; %Range of Fourier components to consider (Vector)
fudge1=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
fudge2=exp(-pi*((mvect/tau_rep-f)/fmax).^2);

dT1=zeros(1,length(mvect))';
dT2=zeros(1,length(mvect))';
boundary=1/sqrt(ws^2+wp^2)*2;
a=-2*pi*boundary;
b=2*pi*boundary;
c=-2*pi*boundary;
d=2*pi*boundary;

for i = 1:length(mvect)
    dT1(i)=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),a,b,c,d,mvect(i)/tau_rep+f);
    dT2(i)=rombint2D(@(ETA,NU,FF) TDTR_3DAni_getTintegral(ETA,NU,FF,k,C,h,wp,Qp,ws),a,b,c,d,mvect(i)/tau_rep-f);
end
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

end

