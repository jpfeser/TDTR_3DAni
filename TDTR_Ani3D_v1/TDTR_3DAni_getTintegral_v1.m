function dT = TDTR_3DAni_getTintegral_v1(fvect,k,C,h,wp,Qp,ws,xoffset,yoffset,xi_min,xi_max,eta_min,eta_max,nnodes);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<14
    nnodes = 35;
end
%initialize some stuff
Nfreq=length(fvect);
InnerInt = zeros(Nfreq,nnodes);
Teta =zeros(Nfreq,nnodes);
dT = zeros(Nfreq);

%%First calculate the nodes/weights for integration
[xi_nodes,xi_weights]=TDTR_Ani3D_lgwt_v1(nnodes,xi_min,xi_max);
[eta_nodes,eta_weights]=TDTR_Ani3D_lgwt_v1(nnodes,eta_min,eta_max);

%do the integral
for i = 1:nnodes
    for j = 1:nnodes
        temp = TDTR_3DAni_getTintegrand_v1(xi_nodes(i),eta_nodes(j),fvect,k,C,h,wp,Qp,ws,xoffset,yoffset);
        InnerInt(:,j) = temp;
    end
    Teta(:,i) = InnerInt*eta_weights;    % do the integral along eta, for all frequencies, at a particular xi
end
dT = Teta*xi_weights; %do the integral along xi, for all frequencies

end

