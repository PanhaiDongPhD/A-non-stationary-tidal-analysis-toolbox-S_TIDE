function [St,Ht,Gt,H,G,coef,xout,ju,coefint]=s_tide_m1(xin,IDP1,IDP2,tides,ntides,n1,dt)
%Haidong Pan(Ocean University of China) 2017/07/01
%email:1027128116@qq.com
% s_tide主程序的修改版本
% A modified version of s_tide
%前n1个分潮使用独立点方案，其余分潮使用经典调和分析
%The first n1 tidal constituents use enhanced harmonic analysis while rest constituents still use classical harmonic analysis
% 
%see s_demo.m for examples
%Reference:Wang, D., H. Pan, G. Jin and X. Lv (2020), Seasonal variation of the main tidal constituents in the Bohai Bay, Ocean Science.
%
disp(['********* Spline Interpolation Tides Analysis(STide)**********  '])
disp(['************* written by Haidong Pan(OUC) ********************  '])
disp(['Independent points of mean sea level:',num2str(IDP1)])
disp(['Independent points of amplitude and phase:',num2str(IDP2)])

load t_constituents
for i=1:ntides
    ju(i)=strmatch(upper(tides(i)),const.name);
    fu(i)=const.freq(ju(i));
end
nobs=length(xin);
gd=find(isfinite(xin(1:nobs)));
ngood=length(gd);
fprintf('   Points used: %d of %d\n',ngood,nobs)

S1 = s_calculate_coefficient(IDP1,nobs);
S2 = s_calculate_coefficient(IDP2,nobs);
disp([' Spline Interpolation coefficients calculated  '])
nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point

t=dt*([1:nobs]'-ceil(nobsu/2));  % Time vector for entire time series,
JJ=n1;

tc=zeros(nobs,IDP1+2*IDP2*JJ+(ntides-JJ)*2,'double');
tc(:,1:IDP1)=S1;
for i=1:nobs
   for j=1:JJ
    for k=1:IDP2
      tc(i,IDP1+(j-1)*IDP2+k)=S2(i,k)*cos((2*pi)*t(i)*fu(j));
      tc(i,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(i,k)*sin((2*pi)*t(i)*fu(j));
    end
   end
end
 for i=1:nobs
   for j=1:ntides-JJ
      tc(i,IDP1+2*IDP2*JJ+j)=cos((2*pi)*t(i)*fu(j+n1));
      tc(i,IDP1+2*IDP2*JJ+ntides-JJ+j)=sin((2*pi)*t(i)*fu(j+n1));
   end
end
%coef=tc(gd,:)\xin(gd);
tci=tc';
[coef,coefint]=regress(tci(:,gd)*xin(gd),tci(:,gd)*tc(gd,:),0.05);

disp([' Least squares soln:use regress '])
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

for j=1:ntides-JJ
    Sa=coef(IDP1+2*IDP2*JJ+j);
    Sb=coef(IDP1+2*IDP2*JJ+j+ntides-JJ);
    H(j)=sqrt(Sa^2+Sb^2);
    G(j)=atand(Sb/Sa);
end
G(Sa<0)=G(Sa<0)+180;

xout=tc*coef;
varx=cov(xin(gd));varxp=cov(xout(gd));
disp(['var(x)= ',num2str(varx),'   var(xp)= ',num2str(varxp)]);
disp(['percent var predicted/var original=',num2str(100*varxp/varx),'%']);
disp(['RMSE=',num2str(sqrt(nanmean((xout(gd)-xin(gd)).^2)))]);
disp(['Max error=',num2str(max(abs(xout(gd)-xin(gd))))]);




