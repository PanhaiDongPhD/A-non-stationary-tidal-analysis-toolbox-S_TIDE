function [St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint,aa,bb,sigma_aa,sigma_bb]=s_tide_m3(xin,IDP1,IDP2,tides,ntides,dt,interpolation,method,latitude,stime)
%Haidong Pan(Ocean University of China) 2020/01/22   潘海东（中国海洋大学）
%email:1027128116@qq.com      panhaidong_phd@qq.com     
% S_TIDE users QQ group number (QQ群号）485997351
% a modified version of s_tide.m (for satellite data)
% s_tide.m主程序的修改版本（处理卫星数据）
% s_tide_m3.m includes nodal factors f and u in harmonic analysis
% s_tide_m3.m  在潮汐调和分析模型里考虑了交点因子和订正角f和u
% Also note: All versions of s_tide can be downloaded from my researchgate：
% https://www.researchgate.net/project/A-non-stationary-tidal-analysis-toolbox-S-TIDE
% 
% new input: latitude (the latitude of observations)
%             stime (the start date of observations)     
%S_TIDE v1.16  

disp(['********** Enhanced Harmonic Analysis(EHA)************  '])
disp(['************* written by Haidong Pan(OUC) ********************  '])
disp(['Independent points of instantaneous sea level:',num2str(IDP1)])
disp(['Independent points of amplitude and phase:',num2str(IDP2)])
siz=size(xin);
if siz(1)<2
    xin=xin';
end
load t_constituents
if iscell(tides)==1
for i=1:ntides
    ju(i)=strmatch(upper(tides(i)),const.name);
    fu(i)=const.freq(ju(i));
end
else
for i=1:ntides
   ju(i)=0;
   fu(i)=tides(i);
end
end
nobs=length(xin);
gd=find(isfinite(xin(1:nobs)));
ngood=length(gd);
fprintf('   Points used: %d of %d\n',ngood,nobs)

ntime=nobs;
ff=zeros(ntides,ntime,'double');
uu=ff;
vv=ff;
for i=1:ntime
days=datenum(stime)+(i-1)*dt/24.0;
[v,u,f]=t_vuf('nodal',days,ju,latitude);
ff(:,i)=f(:);
uu(:,i)=360*u;
vv(:,i)=360*v;
end
nobsu=ntime-rem(ntime-1,2);

for i=1:ntime
vv(:,i)=vv(:,ceil(nobsu/2));
%ff(:,i)=ff(:,ceil(nobsu/2));
%uu(:,i)=uu(:,ceil(nobsu/2));
end
ff=ff';
uu=uu'*2*pi/360;
vv=vv'*2*pi/360;



if strcmp(interpolation,'spline')|strcmp(interpolation,'Spline')
S1 = s_calculate_coefficient(IDP1,nobs); %use spline interpolation 
S2 = s_calculate_coefficient(IDP2,nobs);
disp([' Spline Interpolation coefficients (weights) calculated  '])
end
if strcmp(interpolation,'linear')|strcmp(interpolation,'Linear')
S1 = l_calculate_coefficient(IDP1,nobs);  %use linear interpolation
S2 = l_calculate_coefficient(IDP2,nobs);
disp([' Linear Interpolation coefficients (weights) calculated  '])
end
if strcmp(interpolation,'sinc')|strcmp(interpolation,'Sinc')
S1 = sinc_calculate_coefficient(IDP1,nobs); %use sinc interpolation 
S2 = sinc_calculate_coefficient(IDP2,nobs);
disp([' Sinc Interpolation coefficients (weights) calculated  '])
end


nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point

t=dt*([1:nobs]'-ceil(nobsu/2));  % Time vector for entire time series,
JJ=length(fu);


if (nobs*(IDP1+2*IDP2*JJ))<8760*1700*10
 tc=zeros(nobs,IDP1+2*IDP2*JJ,'double');
 tc(:,1:IDP1)=S1;
 %for i=1:nobs
   for j=1:JJ
    for k=1:IDP2
      tc(1:nobs,IDP1+(j-1)*IDP2+k)=S2(1:nobs,k).*ff(1:nobs,j).*cos((2*pi)*t(1:nobs)*fu(j)+uu(1:nobs,j)+vv(1:nobs,j));
      tc(1:nobs,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(1:nobs,k).*ff(1:nobs,j).*sin((2*pi)*t(1:nobs)*fu(j)+uu(1:nobs,j)+vv(1:nobs,j));
    end
   end
 %end
 tci=tc';
 %coef=(tci(:,gd)*tc(gd,:))\(tci(:,gd)*xin(gd));
if strcmp(method,'ols')|strcmp(method,'Ols')|strcmp(method,'OLS')
   [coef]=regress(tci(:,gd)*xin(gd),tci(:,gd)*tc(gd,:));
   disp([' ***Ordinary Least Squares used!*** '])
end
 if strcmp(method,'robustfit')|strcmp(method,'Robustfit')
   [coef]=robustfit(tc(gd,:),xin(gd),'cauchy',2.385,'off');
   disp([' ***Robust fit weighting method used!*** '])
end
 xout=tc*coef;
 resid=xin(gd)-xout(gd);
sig_res=resid'*resid/(ngood-IDP1-2*IDP2*ntides);
sig_coef=sqrt(inv(tci*tc)*sig_res);
se=diag(sig_coef);se=se';
else
       disp([' ***subblock technology used!*** '])
    nsub=round(nobs/2);  %subblock technology from T_TIDE,original matirx was divided into 50 subblocks, if your records are too long, you need to increase subblocks!
     lhs=zeros(IDP1+2*IDP2*JJ,IDP1+2*IDP2*JJ,'double'); rhs=zeros(IDP1+2*IDP2*JJ,1,'double');
    for j1=1:nsub:ngood                       
      j2=min(j1 + nsub - 1,ngood);
      tc=zeros(j2-j1+1,IDP1+2*IDP2*JJ,'double');
      tc(:,1:IDP1)=S1(gd(j1:j2),:);
      %for i=j1:j2
        for j=1:JJ
          for k=1:IDP2
             tc(1:j2-j1+1,IDP1+(j-1)*IDP2+k)=S2(gd(j1:j2),k).*ff(gd(j1:j2),j).*cos((2*pi)*t(gd(j1:j2))*fu(j)+uu(gd(j1:j2),j)+vv(gd(j1:j2),j));
             tc(1:j2-j1+1,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(gd(j1:j2),k).*ff(gd(j1:j2),j).*sin((2*pi)*t(gd(j1:j2))*fu(j)+uu(gd(j1:j2),j)+vv(gd(j1:j2),j));
          end
        end
      %end
      rhs=rhs + tc'*xin(gd(j1:j2));
      lhs=lhs + tc'*tc;

     
    end;
 % coef=lhs\rhs;
[coef]=regress(rhs,lhs);
disp(['******************************'])
disp(['*****Your records are too long! If use robustfit, lots of computation memories and time are needed. Thus,here,we still use OLS****'])
disp(['******************************'])
  for j1=1:nsub:nobs
      j2=min(j1 + nsub - 1,nobs);
       tc=zeros(j2-j1+1,IDP1+2*IDP2*JJ,'double');
      tc(:,1:IDP1)=S1(j1:j2,:);
      %for i=j1:j2
        for j=1:JJ
          for k=1:IDP2
             tc(1:j2-j1+1,IDP1+(j-1)*IDP2+k)=S2(j1:j2,k).*ff(j1:j2,j).*cos((2*pi)*t(j1:j2)*fu(j)+uu(j1:j2,j)+vv(j1:j2,j));
             tc(1:j2-j1+1,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(j1:j2,k).*ff(j1:j2,j).*sin((2*pi)*t(j1:j2)*fu(j)+uu(j1:j2,j)+vv(j1:j2,j));
          end
        end
      %end
      xout(j1:j2)=tc*coef;
  end;
  xout=xout';
   resid=xin(gd)-xout(gd);
   sig_res=resid'*resid/(ngood-IDP1-2*IDP2*ntides);
   sig_coef=sqrt(inv(lhs)*sig_res);
   se=diag(sig_coef);se=se';
  
end

%disp([' Least squares soln:use regress '])
S_IDP=coef(1:IDP1);
St=S1*S_IDP;

for j=1:JJ
     Sa_IDP=coef(IDP1+1+IDP2*(j-1):IDP1+IDP2*j);
     Sb_IDP=coef(IDP1+1+IDP2*(j-1+JJ):IDP1+IDP2*(j+JJ));
     aa(j,:)=S2*Sa_IDP;
     bb(j,:)=S2*Sb_IDP;
end

Ht=sqrt(aa.^2+bb.^2);
Gt=atand(bb./aa);
Gt(aa<0)=Gt(aa<0)+180;
%xout=tc*coef;
varx=cov(xin(gd));varxp=cov(xout(gd));
disp(['var(x)= ',num2str(varx),'   var(xp)= ',num2str(varxp)]);
disp(['percent var predicted/var original=',num2str(100*varxp/varx),'%']);
disp(['RMSE=',num2str(sqrt(nanmean((xout(gd)-xin(gd)).^2)))]);
disp(['Max error=',num2str(max(abs(xout(gd)-xin(gd))))]);

%计算振幅和迟角的95%置信区间  calculate the  95% confidence interval 
for j=1:JJ
     Sa_IDP=se(IDP1+1+IDP2*(j-1):IDP1+IDP2*j);
     Sb_IDP=se(IDP1+1+IDP2*(j-1+JJ):IDP1+IDP2*(j+JJ));
     sigma_aa(j,:)=S2*Sa_IDP';
     sigma_bb(j,:)=S2*Sb_IDP';
end
deri_h=(aa./Ht.*sigma_aa).^2+(bb./Ht.*sigma_bb).^2;
Htint=1.96*sqrt(deri_h);  %assuming student t distribution

deri_g=(bb./Ht./Ht.*sigma_aa).^2+(aa./Ht./Ht.*sigma_bb).^2;
Gtint=1.96*sqrt(deri_g);

S_IDP=se(1:IDP1);
Stint=1.96*(S1*S_IDP');



