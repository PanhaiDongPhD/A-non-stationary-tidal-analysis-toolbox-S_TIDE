function [Hc,Gc,ff,uu,vv]=s_nodal_correction(Ht,Gt,stime,ntime,dt,latitude,ju)
% written by Haidong Pan, 2017-05-02
% Thanks to Daosheng Wang for pointing out the mistake in correcting phase
% 
% nodal correction of the amplitues and phases.
%stime   tidal signal start time
%ntime   total number of the tidal signal (hour)
%dt       sampling interval (hour)

ntides=length(ju);
ff=zeros(ntides,ntime,'double');
uu=ff;
vv=ff;
disp(['******calculating astronomical phase V, the nodal phase modulation U, and the nodal amplitude correction F*********'])
disp(['***********************accurate nodal/satellite corrections need lots of time!*************************************'])

r0=floor(ntime/2);
r1=floor(ntime/10);
r2=floor(ntime/30);
r3=floor(ntime/60);
r4=floor(ntime/90);

for i=1:ntime
  if i==r0|i==r1|i==r2|i==r3|i==r4
     fprintf('   Numbers calculated: %d of %d\n',i,ntime)
  end
days=datenum(stime)+(i-1)*dt/24.0;
[v,u,f]=t_vuf('nodal',days,ju,latitude);
ff(:,i)=f(:);
uu(:,i)=360*u;
vv(:,i)=360*v;
end
nobsu=ntime-rem(ntime-1,2);
for i=1:ntime
vv(:,i)=vv(:,ceil(nobsu/2));
end
Hc=Ht./ff;
Gc=Gt+uu+vv;
