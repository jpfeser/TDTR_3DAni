function [integral]=trap2D(f,X,Y,Z) 
% computes the integratal I(z)=int(f(x,y,z)) over the square
%   a<x<b, c<y<d using the trapezoidal rule
%   X,Y,Z are an nd meshgrid

[Nx,Ny,Nz]=size(X);
h=-1/(Nx-1)*(X(1,1,1)-X(Nx,1,1)); %distance between nodes, x direction
k=-1/(Ny-1)*(Y(1,1,1)-Y(1,Ny,1)); %distance between nodes, y direction

%generate the weights for the trapezoidal sum (in 2D)
weights=ones(Nx,Ny,Nz)*4;
weights(2:Nx-1,1,:)=2;
weights(2:Nx-1,Ny,:)=2;
weights(1,2:Ny-1,:)=2;
weights(Nx,2:Ny-1,:)=2;
weights(1,1,:)=1;
weights(Nx,1,:)=1;
weights(1,Ny,:)=1;
weights(Nx,Ny,:)=1;

weights=h*k/4*weights.*f(X,Y,Z); %save some memory by reusing weights

integral=squeeze(sum(sum(weights))); %sum all the elements to get the integral;
%use squeeze to get rid of non-useful dimensions

end

