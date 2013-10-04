%This line should be run before entering this script
%[ratio_data,dR_data] = TDTR_Ani3D_EvalSweep(tdelay,lambda_tensor,C,h,f,tau_rep,wp,Qp,ws,TCR,Offset_vect)

%initialize some stuff, only allow offset value.
M_tdelay = length(tdelay);
L_C = length(C);

dR_temp = zeros(M_tdelay);
ratio_temp = zeros(M_tdelay);
S_C = zeros(M_tdelay,L_C);
S_L1 = zeros(M_tdelay,L_C);
S_L2 = zeros(M_tdelay,L_C);
S_L3 = zeros(M_tdelay,L_C);
S_h  = zeros(M_tdelay,L_C);
S_r = zeros(M_tdelay);

for ii=1:length(C)
    %-------------Specific heat-----------------
    Ctemp=C;
    Ctemp(ii) = C(ii)*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,Ctemp,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*C(ii);
    Denom = Ctemp(ii) - C(ii);
    S_C(:,ii) = Num./(ratio_temp*Denom);
    %-------------Thermal Conductivity Tensor ---------
    %-----------(kx)----------
    lambdatemp=lambda_tensor;
    lambdatemp(ii,1)=lambda_tensor(ii,1)*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambdatemp,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*lambdatemp(ii,1);
    Denom = lambdatemp(ii,1) - lambda_tensor(ii,1);
    S_L1(:,ii) =  Num./(ratio_temp*Denom);
        %-----------(ky)----------
    lambdatemp=lambda_tensor;
    lambdatemp(ii,2)=lambda_tensor(ii,2)*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambdatemp,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*lambdatemp(ii,2);
    Denom = lambdatemp(ii,2) - lambda_tensor(ii,2);
    S_L2(:,ii) =  Num./(ratio_temp*Denom);
    %-----------(kz)----------
    lambdatemp=lambda_tensor;
    lambdatemp(ii,3)=lambda_tensor(ii,3)*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambdatemp,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*lambdatemp(ii,3);
    Denom = lambdatemp(ii,3) - lambda_tensor(ii,3);
    S_L3(:,ii) =  Num./(ratio_temp*Denom);
    %-------------Layer Thickess---------------
    %-------------Layer Thickess---------------
    htemp=h;
    htemp(ii) = h(ii)*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,htemp,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*h(ii);
    Denom = htemp(ii) - h(ii);
    S_h(:,ii) =  Num./(ratio_temp*Denom);
    %--------------------------------------------
end

%assumes same size pump/probe to save time
    wtemp=wp_x;
    wtemp = wp_x*1.01;
    [Vout_waste,Vin_waste,ratio_temp] = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wtemp,wtemp,Qp,wtemp,wtemp,TCR,xoffset,yoffset);
    delta = ratio_temp-ratio_data;
    Num = delta*wp_x;
    Denom = wtemp - wp_x;
    S_r =  Num./(ratio_temp*Denom);
    
    %% Make plots
        figure
        semilogx(tdelay,S_C,'*')
        hold on
        semilogx(tdelay,S_L1,'s')
        semilogx(tdelay,S_L2,'d')
        semilogx(tdelay,S_L3,'o')
        semilogx(tdelay,S_h,'x')
        semilogx(tdelay,S_r,'-')
        Cplotlab=strcat(' C ',int2str((1:length(C))'));
        L1plotlab=strcat('L1 ',int2str((1:length(C))'));
        L2plotlab=strcat('L2 ',int2str((1:length(C))'));
        L3plotlab=strcat('L3 ',int2str((1:length(C))'));
        hplotlab=strcat(' h ',int2str((1:length(C))'));
        legend([Cplotlab;L1plotlab;L2plotlab;L3plotlab;hplotlab;' R '])
        