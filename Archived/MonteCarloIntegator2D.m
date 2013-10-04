function [integral] = MonteCarloIntegator2D(f,a,b,c,d,z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Nmax=15;
oldintegral=1e27*ones(length(z),1)
count=0;
for n=5:Nmax
    count=count+1;
    Ns=2^n;
    x=(b-a)*rand(Ns,1)+a; %generate points in range
    y=(d-c)*rand(Ns,1)+c;
    [X,Z]=ndgrid(x,z);
    [Y,Z]=ndgrid(y,z);
    integral=(b-a)*(d-c)/(Ns)*squeeze(sum(f(X,Y,Z)));
    if count~=1
        integral = (2*integral+1*oldintegral)/3;
    end
    relerr = max(abs((integral-oldintegral))/abs(oldintegral));
    if relerr<1e-4
        return
    elseif n==Nmax
        fprintf('too many function evaluations\n')
        return
    end
    relerr;
    oldintegral=integral;
end

