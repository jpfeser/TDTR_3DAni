%The main program tries to minimize "Z" by optimizing the variable(s) X
%This program Lets you:
%   1) Define the vector X: example, if lambda(3) is what you want to solve for,
%   then set X=lambda(3)....if you whish to simulatenous solve for more than one
%   variable, you can just define multiple variables (eg. X(1)=lambda(3), X(2)=lambda(4))
%
%   2) Define the fit.  Typically, this is the sum of the squares of the
%   residuals, but you might want to weight these by the sensitivity,
%   particularly if you don't intend to calculate the errorbars!
%[ratio_data,deltaR_data]=TDTR_3DAni_getREFL_wOffset(tdelay,lambda_tensor,C,h,f,tau_rep,wp,Qp,ws,TCR,xoffset,yoffset)
                                                    
function [Z,ratio_model,deltaR_model]=TDTR_3DAni_FIT(X,ratio_data,tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset)
%Define the variables to be fit
%lambda_tensor is of form lamba_tensor(layer_number,direction_index)
%direction index is 1 (kxx), 2 (kyy), 3 (kzz), 4 (kxy), 5 (kxz), 6 (kyz)
lambda_tensor(3,3)=X(1); 
lambda_tensor(2,3)=X(2);

[Vout_model,Vin_model,ratio_model] =TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
%Uncomment the next three lines to see the non-linear optimization in
%action!
figure(10)
semilogx(tdelay,ratio_data,'ob',tdelay,ratio_model,'g'); 
axis([100e-12 10e-9 0 max([ratio_data;ratio_model])])
pause(0.1)
X
res=(ratio_model-ratio_data).^2;
Z=sum(res);