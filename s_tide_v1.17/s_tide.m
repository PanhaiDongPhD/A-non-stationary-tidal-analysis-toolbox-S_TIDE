function [St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(xin,IDP1,IDP2,tides,ntides,dt,interpolation,method)
%Haidong Pan(Ocean University of China) 2017/05/02   潘海东（中国海洋大学）
%email:1027128116@qq.com      panhaidong_phd@qq.com     
% S_TIDE users QQ group number (QQ群号）485997351
%
%the output:
%St   time-varying sub-tidal oscillations          
%Ht   time-varying amplitudes (contain nodal variations)
%Gt   time-varying phases (contain nodal variations)
%xout   tidal prediction, xin-xout is the residual 
%Stint  the 95% confidence interval for St
%Htint  the 95% confidence interval for Ht
%Gtint  the 95% confidence interval for Gt
%
%
%the input:
%xin    the input tidal signal 
%IDP1   the IP number for sub-tidal oscillations 
%IDP2   the IP number for constituents
%tides    the name of used constituents (you can also use nontidal frequency, for examples see s_demo.m)
%ntides   the number of used constituents, if tides and ntides are all inputed as
%          strings 'autoselected', tides will be selected according to the
%          length of the data (see demo.m for examples)
%dt       sampling interval (hours)
%interpolation    the interpoaltion method used in IP scheme, users can
%                 choose 'spline', 'linear', or 'sinc'
%method     the fitting method used in inversion, users can choose 'ols'
%(Ordinary Least Squares) or 'robustfit' (Robust fit weighting method)
%
%see s_demo.m for the examples
%
%Please citing:
%1. Pan, H., X. Lv, Y. Wang, P. Matte, H. Chen, and G. Jin (2018), Exploration of Tidal-Fluvial Interaction in the Columbia River Estuary Using S_TIDE, J. Geophys. Res. Ocean., 123(9), 6598–6619, doi:10.1029/2018JC014146.
%2.	Jin, G., H. Pan, Q. Zhang, X. Lv, W. Zhao, and Y. Gao (2018), Determination of Harmonic Parameters with Temporal Variations: An Enhanced Harmonic Analysis Algorithm and Application to Internal Tidal Currents in the South China Sea, J. Atmos. Ocean. Technol., 35(7), 1375-1398
%3. Pawlowicz, R., B. Beardsley, and S. Lentz (2002), Classical tidal harmonic analysis including error estimates in MATLAB using T_TIDE, Comput. Geosci., 28(8), 929–937
%
%Haidong Pan  2017/05/02 - S_TIDE was developed from the widely used T_TIDE to
%                          realize EHA,
%             2018/08/29  - After many corrections, S_TIDE v1.0 matlab toolbox  was published
%             2018/09/08  - S_TIDE v1.1 matlab toolbox  was published,using
%             subblock technology from T_TIDE, S_TIDE needs less
%             computation memory and time than before.
%             2018/09/20  - S_TIDE v1.11 toolbox was published,
%             s_tide_error function was embedded into s_tide function,also
%             the old error estimation algorithm (see Pan et al.2018 for
%             details) was repalced by the new error estimation algorithm
%             (see Foreman and Henry,1989 for details)
%             2018/11/10  -S_TIDE v1.12 toolbox was published, a new
%             interpolation method--the Sinc interpolation was added into
%             S_TIDE.
%             2018/12/16 -S_TIDE v1.13 toolbox was published,a new function
%             s_dlrm.m was added. s_dlrm.m uses Dynamic Linear Regression Model(DLRM), 
%             which can calculate time-varying regression coefficients (see
%             Pan, 2018 for details). See s_demo.m for examples.
%             2019/01/12 -S_TIDE v1.14 was published, robustfit function
%             was added, which can reduce the contribution of high-leverage data points, in order
%             to improve the overall fit.  
%             2019/06/25  -S_TIDE v1.15 was published, non-tidal frequency
%             can be added by users. Also, when IDP=1, S_TIDE outputs
%             constant harmonic parameters and MWL.When IDP=2,S_TIDE calculate the
%             linear trend of harmonic parameters and MWL.
%             2020/01/13  - S_TIDE v1.16 was published, a modified version
%             of s_tide.m was added, namely (s_tide_m1.m). S_tide_m1.m assumes that the first n1 tidal constituents use enhanced harmonic analysis 
%             while rest constituents still use classical harmonic analysis
%             2020/05/10-S_TIDE v1.17 was published, when IDP2=1, S_TIDE
%             can automatically select tides according to the length of data
% 
% Citing:
% Foreman and Henry (1989)  The harmonic analysis of tidal model time series
% Pan et al.(2018)  Exploration of tidal-fluvial interaction in the Columbia River estuary using S_TIDE
% Pan, H.(2018) Application Examples of S_TIDE Toolbox. Technical Report 2018-11. Key Laboratory of Physical Oceanography, Ocean University of China, Qingdao, China, 20pp. https://www.researchgate.net/publication/329091915_Application_Examples_of_S_TIDE_Toolbox
% 
% Also note: All versions of s_tide can be downloaded from my researchgate：
% https://www.researchgate.net/project/A-non-stationary-tidal-analysis-toolbox-S-TIDE
%
% S_TIDE v1.17

disp(['********** Enhanced Harmonic Analysis(EHA)************  '])
disp(['************* written by Haidong Pan(OUC) ********************  '])
disp(['Independent points of mean sea level:',num2str(IDP1)])
disp(['Independent points of amplitude and phase:',num2str(IDP2)])
siz=size(xin);
if siz(1)<2
    xin=xin';
end
nobs=length(xin);
gd=find(isfinite(xin(1:nobs)));
ngood=length(gd);
fprintf('   Points used: %d of %d\n',ngood,nobs)
nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point
t=dt*([1:nobs]'-ceil(nobsu/2));  % Time vector for entire time series,

load t_constituents
if iscell(tides)==1
for i=1:ntides
    ju(i)=strmatch(upper(tides(i)),const.name);
    fu(i)=const.freq(ju(i));
end

elseif isstr(tides)==1
[nameu,fu,ju]=constituents(1/(dt*nobsu));   
ntides=length(nameu);
elseif isnumberic(tides)==1
for i=1:ntides
   ju(i)=0;
   fu(i)=tides(i);
end
end

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

JJ=length(fu);

if (nobs*(IDP1+2*IDP2*JJ))<8760*300
 tc=zeros(nobs,IDP1+2*IDP2*JJ,'double');
 tc(:,1:IDP1)=S1;
 %for i=1:nobs
   for j=1:JJ
    for k=1:IDP2
      tc(1:nobs,IDP1+(j-1)*IDP2+k)=S2(1:nobs,k).*cos((2*pi)*t(1:nobs)*fu(j));
      tc(1:nobs,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(1:nobs,k).*sin((2*pi)*t(1:nobs)*fu(j));
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
   [coef,status]=robustfit(tc(gd,:),xin(gd),'cauchy',2.385,'off');
   disp([' ***Robust fit weighting method used!*** '])
end
 xout=tc*coef;
 resid=xin(gd)-xout(gd);
sig_res=resid'*resid/(ngood-IDP1-2*IDP2*ntides);
sig_coef=sqrt(inv(tci*tc)*sig_res);
se=diag(sig_coef);se=se'; % correlated bivariate white-noise model

if isstr(tides)==1 %uncorrelated bivariate coloured-noise model
  xres=xin-xout; % and the residuals!
  xr=fixgaps(xres);
  [NP,NM]=noise_realizations(xr(isfinite(xr)),fu,dt,300);
  se(IDP1+1:IDP1+JJ)=std(real(2*NM(:,2:end)),0,2);
  se(IDP1+1+JJ:IDP1+2*JJ)=std(imag(2*NM(:,2:end)),0,2);
end
else
    nsub=round(nobs/50);  %subblock technology from T_TIDE,original matirx was divided into 50 subblocks, if your records are too long, you need to increase subblocks!
     lhs=zeros(IDP1+2*IDP2*JJ,IDP1+2*IDP2*JJ,'double'); rhs=zeros(IDP1+2*IDP2*JJ,1,'double');
    for j1=1:nsub:ngood                       
      j2=min(j1 + nsub - 1,ngood);
      tc=zeros(j2-j1+1,IDP1+2*IDP2*JJ,'double');
      tc(:,1:IDP1)=S1(gd(j1:j2),:);
      %for i=j1:j2
        for j=1:JJ
          for k=1:IDP2
             tc(1:j2-j1+1,IDP1+(j-1)*IDP2+k)=S2(gd(j1:j2),k).*cos((2*pi)*t(gd(j1:j2))*fu(j));
             tc(1:j2-j1+1,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(gd(j1:j2),k).*sin((2*pi)*t(gd(j1:j2))*fu(j));
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
             tc(1:j2-j1+1,IDP1+(j-1)*IDP2+k)=S2(j1:j2,k).*cos((2*pi)*t(j1:j2)*fu(j));
             tc(1:j2-j1+1,IDP1+IDP2*JJ+(j-1)*IDP2+k)=S2(j1:j2,k).*sin((2*pi)*t(j1:j2)*fu(j));
          end
        end
      %end
      xout(j1:j2)=tc*coef;
  end;
  xout=xout';
   resid=xin(gd)-xout(gd);
   sig_res=resid'*resid/(ngood-IDP1-2*IDP2*ntides-1);
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
Gtint=1.96*sqrt(deri_g)*180/pi;

S_IDP=se(1:IDP1);
Stint=1.96*(S1*S_IDP');
%-----------------Output results to the screen---------------------------------------
 if isstr(tides)==1
    fprintf('\n     tidal amplitude and phase with 95%% CI estimates\n');
    fprintf('\ntide   freq       amp     amp_err    pha    pha_err     snr\n');
    snr=(Ht(:,1)./Htint(:,1)).^2;
    for k=1:length(fu)
      if snr(k)>=2, fprintf('*'); else, fprintf(' '); end
      fprintf('%s %9.7f %9.4f %8.3f %8.2f %8.2f %8.2g\n',nameu(k,:),fu(k),Ht(k,1),Htint(k,1),Gt(k,1),Gtint(k,1),snr(k));
    end
 end

function [nameu,fu,ju]=constituents(minres)
% modified from T_TIDE
%2020/05/10
%if minres>1/(18.6*365.25*24)                        % Choose only resolveable pairs for short
 load t_constituents
 ju=find(const.df>=minres);
disp(['   number of standard constituents used: ',int2str(length(ju))])
nameu=const.name(ju,:);
fu=const.freq(ju);
% Check if neighboring chosen constituents violate Rayleigh criteria.
jck=find(diff(fu)<minres);
if (~isempty(jck))
   disp(['  Warning! Following constituent pairs violate Rayleigh criterion']);
   for ick=1:length(jck)
   disp(['     ',nameu(jck(ick),:),'  ',nameu(jck(ick)+1,:)]);
   end
end

function [NP,NM]=noise_realizations(xres,fu,dt,nreal)
% modified from T_TIDE
  [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt);
  Pxcave=zeros(size(Pxcave));  %% For comparison with other technique!
  %fprintf('**** Assuming no covariance between u and v errors!*******\n');
nfband=size(fband,1);

Mat=zeros(4,4,nfband);
for k=1:nfband

  % The B matrix represents the covariance matrix for the vector
  % [Re{ap} Im{ap} Re{am} Im{am}]' where Re{} and Im{} are real and
  % imaginary parts, and ap/m represent the complex constituent 
  % amplitudes for positive and negative frequencies when the input
  % is bivariate white noise. For a flat residual spectrum this works 
  % fine.
 
  % This is adapted here for "locally white" conditions, but I'm still
  % not sure how to handle a complex sxy, so this is set to zero
  % right now.
  
  p=(Pxrave(k)+Pxiave(k))/2;
  d=(Pxrave(k)-Pxiave(k))/2;
  sxy=Pxcave(k);
  
  B=[p    0   d   sxy;
     0    p  sxy  -d;
     d   sxy  p    0
     sxy -d   0    p];

  % Compute the transformation matrix that takes uncorrelated white 
  % noise and makes noise with the same statistical structure as the 
  % Fourier transformed noise.
  [V,D]=eig(B);
  Mat(:,:,k)=V*diag(sqrt(diag(D)));
end

% Generate realizations for the different analyzed constituents.

N=zeros(4,nreal);
NM=zeros(length(fu),nreal);
NP=NM;
for k=1:length(fu)
  l=find(fu(k)>fband(:,1) & fu(k)<fband(:,2));
  N=[zeros(4,1),Mat(:,:,l)*randn(4,nreal-1)];
  NP(k,:)=N(1,:)+i*N(2,:);
  NM(k,:)=N(3,:)+i*N(4,:);
end


function [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt)
% RESIDUAL_SPECTRUM: Computes statistics from an input spectrum over
% a number of bands, returning the band limits and the estimates for
% power spectra for real and imaginary parts and the cross-spectrum.          
%
% Mean values of the noise spectrum are computed for the following 
% 8 frequency bands defined by their center frequency and band width:
% M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd; 
% M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd). 
%from T_TIDE
% S. Lentz  10/28/99
% R. Pawlowicz 11/1/00
% Version 1.0

% Define frequency bands for spectral averaging.
fband =[.00010 .00417;
        .03192 .04859;
        .07218 .08884;
        .11243 .12910;
        .15269 .16936;
        .19295 .20961;
        .23320 .25100;
        .26000 .29000;
        .30000 .50000];

% If we have a sampling interval> 1 hour, we might have to get
% rid of some bins.
%fband(fband(:,1)>1/(2*dt),:)=[];

nfband=size(fband,1);
nx=length(xres);

% Spectral estimate (takes real time series only).


% Matlab has changed their spectral estimator functions
% To match the old code, I have to divide by 2*dt. This is because
% 
%  PSD*dt  is two-sided spectrum in units of power per hertz.
%
%  PWELCH is the one-sided spectrum in power per hertz
%
%  So PWELCH/2 = PSD*dt


%[Pxr,fx]=psd(real(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
%[Pxi,fx]=psd(imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
%[Pxc,fx]=csd(real(xres),imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.


[Pxr,fx]=pwelch(real(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
Pxr=Pxr/2/dt;
[Pxi,fx]=pwelch(imag(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxi=Pxi/2/dt;
[Pxc,fx]=cpsd(real(xres),imag(xres),[],[],nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxc=Pxc/2/dt;

df=fx(3)-fx(2);
Pxr(round(fu./df)+1)=NaN ; % Sets Px=NaN in bins close to analyzed frequencies
Pxi(round(fu./df)+1)=NaN ; % (to prevent leakage problems?).
Pxc(round(fu./df)+1)=NaN ; 

Pxrave=zeros(nfband,1);
Pxiave=zeros(nfband,1);
Pxcave=zeros(nfband,1);
% Loop downwards in frequency through bands (cures short time series
% problem with no data in lowest band).
%
% Divide by nx to get power per frequency bin, and multiply by 2
% to account for positive and negative frequencies.
%
for k=nfband:-1:1
   jband=find(fx>=fband(k,1) & fx<=fband(k,2) & isfinite(Pxr));
   if any(jband)
     Pxrave(k)=mean(Pxr(jband))*2/nx;
     Pxiave(k)=mean(Pxi(jband))*2/nx;
     Pxcave(k)=mean(Pxc(jband))*2/nx;
   elseif k<nfband
     Pxrave(k)=Pxrave(k+1);   % Low frequency bin might not have any points...
     Pxiave(k)=Pxiave(k+1);   
     Pxcave(k)=Pxcave(k+1);   
   end
end
function y=fixgaps(x)
% FIXGAPS: Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN in the input time 
% series (may be complex), but ignores trailing and leading NaNs.
%from T_TIDE
% R. Pawlowicz 11/6/99
% Version 1.0
y=x;
bd=isnan(x);
gd=find(~bd);
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
y(bd)=interp1(gd,x(gd),find(bd)); 


