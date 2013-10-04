%These first 10 lines of code are not important
keepdir=0;
if keepdir==1
    save('bla.mat','filename','pathname','keepdir')
else
    save('bla.mat','keepdir')
end
clear all
pause(0.1)
load('bla.mat','filename','pathname','keepdir');
tic %Start the timer

%-------------TYPE THERMAL SYSTEM PARAMTERS HERE--------------
abslayer =17;
lambda=[17 0.127 80]; %W/m-K
phi = ones(size(lambda)); %radial anisotropic ratio, usually = 1 for isotropic.
lambda_tensor = MakeRadial_ktensor(lambda,phi); %initialize tensor assuming radial symmetry in all layers
%lambda_tensor is of form lamba_tensor(layer_number,direction_index)
%direction index is 1 (kxx), 2 (kyy), 3 (kzz), 4 (kxy), 5 (kxz), 6 (kyz)
lambda_tensor(3,1)=140;
C=[2.65 0.1 3.75]*1e6; %J/m^3-K
h =[65 1 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)

%r=1.0e-6;
wp_x=25e-6; %pump 1/e^2 radius, m
ws_x=25e-6; %probe 1/e^2 radius, m
wp_y = wp_x;
ws_y = ws_x;

f=1.1e6; %laser Modulation frequency, Hz

xoffset=0;
yoffset=0;
Offset = [xoffset,yoffset]; %concentric

tau_rep=1/80e6; %laser repetition period, s
Qp=30e-3; %laser power (Watts) . . . only used for amplitude est.
TCR=1e-4; %coefficient of thermal reflectance . . . only used for amplitude est.
%------------- THERMAL SYSTEM PARAMTERS END HERE--------------

%choose time delays for Sensitivity plots
tdelay=logspace(log10(500e-12),log10(4e-9),20)'; %vector of time delays (used to generate sensitivity plots)

%Choose range of time delays to fit, sec
tdelay_min=500e-12;
tdelay_max=4e-9;

%----------------------------------PROGRAM OPTIONS BEGIN--------
%Generate Sensitivity Plots?
senseplot=1;
%Import Data? 0 for no, 1 for yes
importdata=0;
%If so, which variable(s) are you fitting for?
Xguess=[lambda_tensor(3,3),lambda_tensor(2,3)];%, lambda(3)]; %initial guess for solution, could be for simulated data (if nothing imported) or real data
%Calculate Errorbars? 0 for no, 1 for yes (takes longer)
ebar=0;
%----------------------------------PROGRAM OPTIONS END--------

%-----------------Make sensitivity plots--------------
condition1 = (senseplot==1);
condition2 = (senseplot==0)&(ebar==1);
if (condition1||condition2)
    [Vout_data,Vin_data,ratio_data] = TDTR_3DAni_TDTR_Sig(tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset);
    deltaR_data = Vin_data + sqrt(-1)*Vout_data
end

%temporary
%[tdelay([1 end]),ratio_data([1 end]),real(deltaR_data([1 end])),imag(deltaR_data([1 end]))]
%return
%end temporary

if senseplot==1
    TDTR_Ani3D_GetSensitivities %works..verified
end

%% --------------Import Data---------------------------
if importdata==1
    if keepdir==1
        filename=input('enter file name string:\n')
    else
        [filename,pathname]=uigetfile('*.*','Choose data file');
    end
    DATAMATRIX=dlmread(strcat(pathname,filename),'\t');
    tdelay_raw=DATAMATRIX(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
    Vin_raw=DATAMATRIX(:,3); %Units (uV ?)
    Vout_raw=DATAMATRIX(:,4);
    ratio_raw=DATAMATRIX(:,5);
    [tdelay_data,Vin_data] = extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
    [tdelay_data,Vout_data] = extract_interior(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
    [tdelay_data,ratio_data] = extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
    tdelay = tdelay_data;
    
%--------------Perform Fit (skips if no data import)--------------------------
    Xsol=fminsearch(@(X) TDTR_3DAni_FIT(X,ratio_data,tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset),Xguess)
    Xguess=Xsol;
    tdelay=tdelay_data;
    fprintf('Data fit completed\n')
else
    [tdelay_data,Vin_data] = extract_interior(tdelay,real(deltaR_data),tdelay_min,tdelay_max);
    [tdelay_data,Vout_data] = extract_interior(tdelay,imag(deltaR_data),tdelay_min,tdelay_max);
    [tdelay_data,ratio_data] = extract_interior(tdelay,ratio_data,tdelay_min,tdelay_max);
    Xsol = Xguess;
    tdelay = tdelay_data;
end

fprintf('Fitting Solution:\n')
Xsol
toc

%% Save Data
[Z,ratio_model]=TDTR_3DAni_FIT(Xsol,ratio_data,tdelay,lambda_tensor,C,h,f,tau_rep,wp_x,wp_y,Qp,ws_x,ws_y,TCR,xoffset,yoffset)
dlmwrite('last_fit.txt',[tdelay_data*1e12,ratio_model],'delimiter','\t')
%----------------------------------------------------
saveint=input('Want to save results?\n(0=no, 1=yes)\n');
if saveint==1
    save(strcat(pathname,filename(1:end-4),'Results.mat'))
    fignum=203;
    figure(fignum)
    semilogx(tdelay_data,ratio_data,'ko',tdelay_data,ratio_model,'k-')
    xlabel('time delay (s)','FontSize',18)
    ylabel('-Vin/Vout','FontSize',18)
    title('Data Fit')
    legend('experiment','model')
    set(gca,'FontSize',18)
    axis([min(tdelay_data) max(tdelay_data) 0 1.2*max([ratio_data;ratio_model])])
    
    print(fignum,'-depsc',strcat(pathname,filename(1:end-4),'FIT.eps'))
end
%----------------------------------------------------
fprintf('Program Completed\n')
beep
beep
beep






