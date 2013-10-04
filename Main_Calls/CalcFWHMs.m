%% First calculate concentric case
xoffset = 0;
yoffset = 0;

Vout_max = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);

%% Now find FWHM of xoffset direction
yoffset =0;
Lpx = sqrt(lambda_tensor(3,1)/C(3)*1/(2*pi*f));
xmin_guess = sqrt(wp_x^2+Lpx^2)*0.8;
xmax_guess = sqrt(wp_x^2+Lpx^2)*1.3;
x_FWHM(i) = fzero(@(x) 0.5*Vout_max - TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,x,yoffset),[xmin_guess,xmax_guess],optimset('Display','iter','TolX',1e-9))

%% Now find FWHM of yoffset direction
xoffset =0;
Lpy = sqrt(lambda_tensor(3,2)/C(3)*1/(2*pi*f));
ymin_guess = sqrt(wp_y^2+Lpy^2)*0.8;
ymax_guess = sqrt(wp_y^2+Lpy^2)*1.3;
y_FWHM(i) = fzero(@(x) 0.5*Vout_max - TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,x),[ymin_guess,ymax_guess],optimset('Display','iter','TolX',1e-9))

FWHM_ratio(i) = x_FWHM(i)/y_FWHM(i);
eta(i) =  lambda_tensor(3,1)/ lambda_tensor(3,2);