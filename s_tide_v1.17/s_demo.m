% S_DEMO - demonstration of capabilities.
% Short example of capabilities of S_TIDE.
%%
% demonstration of s_calculate_coefficient.m
x=[1 2 3 4 5];y=[0 1 0 -1 0];
xi=1:0.05:5;
IDP1=5;
nobs=length(xi);
S1 = s_calculate_coefficient(IDP1,nobs); %use spline interpolation 
S2 = l_calculate_coefficient(IDP1,nobs);  %use linear interpolation
S3 = sinc_calculate_coefficient(IDP1,nobs);  %use sinc interpolation
yi_spline=S1*y';
yi_linear=S2*y';
yi_sinc=S3*y';
plot(x,y,'ro','MarkerSize',8.0);hold on;plot(xi,yi_spline,'b','linewidth',1.1);
plot(xi,yi_linear,'m','linewidth',1.1);plot(xi,yi_sinc,'k','linewidth',1.1)
set(gca,'Fontsize',15,'linewidth',1.2);
legend('Given','Spline','Linear','Sinc')
ylim([-1.2 1.2])

%%
% demonstration of s_tide.m
load kushiro.mat   % hourly elevation data at Kushiro (Japan) from 1993-01-01 to 2012-12-31, obtained from PSMSL
load imf.mat      % the modes obtained by Empirical Mode Decomposition (EMD) of the kushiro elevation data.
% If IDP is equal to 2, s_tide will calculate the linear trend of water level!
%If IDP is equal to 1, s_tide will calcualte a constant!
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,1,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St);hold on
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,2,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'r')
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,3,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'k')
legend('IP=1','IP=2','IP=3')
%calculate given frequencies(non-tidal)
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,2,21,[1/36,1/72],2,1,'spline','ols');%period(36 hours and 72 hours)
plot(Ht(1,:),'r'); hold on;plot(Ht(2,:),'k')

%
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,21,21,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'k','linewidth',1.2)
hold on
plot(imf(16,:)+imf(15,:)+imf(14,:)+imf(13,:)+imf(12,:),'r','linewidth',1.2)
legend('S\_TIDE','EMD')
%legend('boxoff')
set(gca,'Fontsize',15,'linewidth',1.2)
set(gca,'XTick',0:8760*2:8760*20)
set(gca,'XTickLabel',{'1993','1995','1997','1999','2001','2003','2005','2007','2009','2011','2013'})
ylabel('Sea Level/mm','Fontsize',15)
xlabel('Year','Fontsize',15)
grid on
xlim([0 8760*20])
%plot the amplitude and phase and their 95% confidence intervals
for i=1:2
figure()
plot(Ht(i,:));hold on
plot(Ht(i,:)+Htint(i,:),'r')
plot(Ht(i,:)-Htint(i,:),'r')
end
for i=1:2
figure()
plot(Gt(i,:));hold on
plot(Gt(i,:)+Gtint(i,:),'r')  
plot(Gt(i,:)-Gtint(i,:),'r')
end
%%
% demonstration of s_nodal_correction.m
% You must have installed T_TIDE to run the code below
[NAME,FREQ,TIDECON,XOUT]=t_tide(kushiro(1:8760*3),'interval',1,'rayleigh',['M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2']); %without nodal correction
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*3),4,4,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','ols');
subplot(2,1,1);plot(Ht(1,:),'r');hold on;plot(TIDECON(6,1)*ones(1,8760*3),'k');ylabel('Amplitude'); title('nodal uncorrected')
subplot(2,1,2);plot(Gt(1,:),'r');hold on;plot(TIDECON(6,3)*ones(1,8760*3),'k');ylabel('Phase')
legend('S\_TIDE','T\_TIDE')

[NAME,FREQ,TIDECON,XOUT]=t_tide(kushiro(1:8760*3),'interval',1,'rayleigh',['M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'],'latitude',42.5,'start time',[1993,01,01,00]);
[Hc,Gc,ff,uu,vv]=s_nodal_correction(Ht,Gt,[1993,01,01],8760*3,1,42.5,ju); figure()
subplot(2,1,1);plot(Hc(1,:),'r');hold on;plot(TIDECON(6,1)*ones(1,8760*3),'k');ylabel('M2 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(1,:),'r');hold on;plot(TIDECON(6,3)*ones(1,8760*3),'k');ylabel('M2 Phase')
legend('S\_TIDE','T\_TIDE')

figure()
subplot(2,1,1);plot(Hc(3,:),'r');hold on;plot(TIDECON(4,1)*ones(1,8760*3),'k');ylabel('K1 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(3,:),'r');hold on;plot(TIDECON(4,3)*ones(1,8760*3),'k');ylabel('K1 Phase')
legend('S\_TIDE','T\_TIDE')

figure()
subplot(2,1,1);plot(Hc(7,:),'r');hold on;plot(TIDECON(3,1)*ones(1,8760*3),'k');ylabel('P1 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(7,:),'r');hold on;plot(TIDECON(3,3)*ones(1,8760*3),'k');ylabel('P1 Phase')
legend('S\_TIDE','T\_TIDE')
%%
%extract 18.61 year cycle from monthly K1 tidal amplitudes using S_TIDE
load K1_amp_hilo.mat    
% the IP number for 18.61 year cycle is one which means the amplitude and phase of 18.61 year cycle are constant
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(K1_amp,2,1,[1/(18.61*8760)],1,720,'spline','ols');
figure();plot(K1_amp,'r','linewidth',1.2);hold on;plot(xout,'k','linewidth',1.2)
 mean(Ht)/mean(St)*100  % theoretical value is 11.5%, but the actual value is 11.72%

% If the IP number for 18.61 year cycle is zero, s_tide will only calculate the linear trend
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(K1_amp,2,0,[1/(18.61*8760)],1,720,'spline','ols');
figure();plot(K1_amp,'r','linewidth',1.2);hold on;plot(xout,'k','linewidth',1.2)
%%
% OLS vs robustfit, this example was first added in S_TIDE v1.14
load kushiro.mat   % hourly elevation data at Kushiro (Japan) from 1993-01-01 to 2012-12-31, obtained from PSMSL
load imf.mat 
tic
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*9),10,10,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','ols');
toc
%ols needs 0.874286s
tic
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*9),10,10,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit');
toc
%robustfit needs 20.675073s
% Ols is much faster than robustfit, while their results are nearly same!

%%
%demonstration of s_tide_m1.m
%this example was first added in S_TIDE v1.16
%Reference:Wang, D., H. Pan, G. Jin and X. Lv (2020), Seasonal variation of the main tidal constituents in the Bohai Bay, Ocean Science.
load kushiro.mat 
[St,Ht,Gt,H,G,coef,xout,ju,coefint]=s_tide_m1(kushiro,200,50,{'M2';'K1';'O1';'S2';'N2';'K2';'P1';'Q1'},8,2,1);
plot(Ht(1,:),'k');hold on;plot(Ht(2,:),'r')
legend('M2','K1');xlabel('Time(hour)');ylabel('Amplitude(mm)')
%我们可以清楚的看到M2和K1分潮的年变化和18.61年变化
%We can clearly see the annual and nodal variations in M2 and K1 amplitudes

%%
%demonstration of s_tide_m2.m (the function s_tide_m1.m can be seen as a special case of s_tide_m2.m)
%this example was first added in S_TIDE v1.16
load kushiro.mat 
[St,Ht,Gt,coef,xout,ju,coefint]=s_tide_m2(kushiro,20,50,5,{'M2';'K1';'O1';'S2';'N2';'K2';'P1';'Q1'},8,2,1);
figure();plot(Ht(1,:),'k');hold on;plot(Ht(2,:),'r');plot(Ht(3,:),'b')
legend('M2','K1','O1');xlabel('Time(hour)');ylabel('Amplitude(mm)')
%我们可以清楚的看到M2和K1分潮的年变化和18.61年变化
%We can clearly see the annual and nodal variations in M2 and K1 amplitudes

%%
%demonstration of s_tide_m3.m
%this example was first added in S_TIDE v1.16
load('satellite.mat');
 cons={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'SA';'SSA'};lat=14.1633;
 stime=[1992,10,31,19,43,37.71];
 [St1,Ht1,Gt1,coef,xout1,ju,Stint1,Htint1,Gtint1,aa1,bb1,sigma_aa1,sigma_bb1]=s_tide_m3(satellite,2,1,cons,10,9.9156*24,'spline','ols',lat,stime);
  
 plot(satellite,'*-','linewidth',1.1);hold on;plot(xout1,'r','linewidth',1.1);
 set(gca,'Fontsize',12,'linewidth',1.2)
 xlim([1 101]);set(gca,'XTick',[1 51 101])
 ylabel('Sea Level(m)');xlabel('Time')
  set(gca,'XTickLabel',{'1992/10/31','1994/03/10', '1995/07/19'})
    legend('Observation','Hindcast')
%%
% this example was first added in S_TIDE v1.17
load kushiro.mat %data start from 1993/01/01
[nameu,fu,tidecon,xout]=t_tide(kushiro(1:8767)); % without nodal corrections
[nameu,fu,tidecon,xout]=t_tide(kushiro(1:8767),'latitude',42.5,'start time',[1993,01,01,00]);%with nodal corrections

% S_TIDE results without nodal corrections
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8767),1,1,'autoselected','autoselected',1,'spline','ols');
%compared to T_TIDE,the results are nearly same
nobs=length(kushiro(1:8767));nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point
ctime=datenum(1993,1,1)+(ceil(nobsu/2)-1)/24.0; % using center time for nodal correction
[v,u,f]=t_vuf('nodal',ctime,ju,42.5);
u=360*u;%convert to degree
v=360*v;
Hc=Ht(:,1)./f; %nodal corrected amplitudes
Gc=Gt(:,1)+u+v; %nodal corrected phases


