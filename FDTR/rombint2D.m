function [sol,nfinal]=rombint2D(f,a,b,c,d,z) 
%rombint2D(f,a,b,c,d,lengthf) integrated the fuction f(x,y,z) over the square
%   a<x<b, c<y<d, and may carry an extra (unintegrated) dimension z or not,
%   so that f=f(x,y,z) or f(x,y)

NN=length(z);
relerr=1;
abserr=1;
limit=1e-6;
nmax=25; %max power for number of points in each grid direction 2^25 = 33 million points
mmax=25;


%integrate over x,y
R=zeros(nmax+1,mmax+1,NN);
n=0;
m=0;
[X,Y,Z]=ndgrid([a b],[c d],z);

R(1:NN,n+1,m+1)=0.25*(b-a)*(d-c)*(sum(sum(f(X,Y,Z)))); %first evaluation

while max(relerr)>limit
    n=n+1;
    % do trap rule on next set of points;
    x=linspace(a,b,2^n+1);
    y=linspace(c,d,2^n+1);
    [X,Y,Z]=ndgrid(x,y,z);
    %
    R(n+1,1,1:NN)=trap2D(f,X,Y,Z);
    %now do romberg extrapolation
    for m=1:n
       fourm=4^m;
       R(n+1,m+1,1:NN)=1/(fourm-1)*(fourm*R(n+1,m,1:NN)-R(n,m,1:NN)); 
    end
    R(n+1,m+1,1:NN);
    relerr=abs(1-R(n,m,1:NN)./R(n+1,m+1,1:NN));
    abserr=abs(R(n+1,m+1,1:NN)-R(n,m,1:NN));
    if n==nmax
        fprintf('max iterations reached')
        sol=R(n+1,m+1,1:NN);
        nfinal=n;
        return
    end
end
sol=R(n+1,m+1,1:NN);
nfinal=n;
end

