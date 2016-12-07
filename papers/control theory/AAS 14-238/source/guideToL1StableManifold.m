function guideToL1StableManifold 
close all; clear; clc

global mue acmd C;
mue=0.01215; C=3.16; bon=1-mue;
options = odeset('RelTol',1e-11,'AbsTol',1e-8);
figure(1)
% re = .15;
% rm = .05;
% [xe ye ze] = sphere;
hold on;
% surf(re*xe-mue,re*ye,re*ze)
% surf(rm*xe+(1-mue),rm*ye,rm*ze)
primaryBodyPos = [-mue;0;0];
secondaryBodyPos = [1-mue;0;0];
plot3(primaryBodyPos(1),primaryBodyPos(2),primaryBodyPos(3),'kx');
plot3(secondaryBodyPos(1),secondaryBodyPos(2),secondaryBodyPos(3),'kx');

xlabel('Synodic-x');
ylabel('Synodic-y');
zlabel('Synodic-z');

% Equilibrium Points
y1=@(x)x-(1-mue)/(x+mue)^2+mue/(x-1+mue)^2;Lp1=fzero(y1,0);
y2=@(x)x-(1-mue)/(x+mue)^2-mue/(x-1+mue)^2;Lp2=fzero(y2,0);
y3=@(x)x+(1-mue)/(x+mue)^2+mue/(x-1+mue)^2;Lp3=fzero(y3,0);
Lp4x=0.5-mue; Lp4y=0.5*sqrt(3);
Lp5x=0.5-mue; Lp5y=-0.5*sqrt(3);
plot(Lp1,0,'r.')
plot(Lp2,0,'r.')
plot(Lp3,0,'r.')
plot(Lp4x,Lp4y,'g.')
plot(Lp5x,Lp5y,'g.')
%--------------------------------------------------------------
% Halo around L1
x10=0.869581; x30=0; x50=-0.05;
x20=0; x40=-Velo(C,x10,x30,x50); x60=0;
[t,x]=ode45(@cr3bp_3b, [0 3], [x10 x20 x30 x40 x50 x60],options);
hx=x(:,1); hy=x(:,3); hz=x(:,5);
hvx=x(:,2); hvy=x(:,4); hvz=x(:,6);n=4;

%plot halo orbit
plot3(hx,hy,hz,'k') %L1
manfldInd = 300;
for i=1:53
    x10=hx(n*i)-0.0001; x30=hy(n*i); x50=hz(n*i);
    x20=hvx(n*i); x40=hvy(n*i); x60=hvz(n*i);
    [t,x]=ode45(@cr3bp_3b, [0 -6], [x10 x20 x30 x40 x50 x60 ],options);m=length(t);
    for ii=1:m
        if x(ii,1)<0
            if x(ii,3)>0
                son=ii;
                break
            end
        end
    end
    if (i == 10) %arbitrary point on manifold
        xf0 = x(manfldInd,1);
        vfx0 = x(manfldInd,2);
        yf0 = x(manfldInd,3);
        vfy0 = x(manfldInd,4);
        zf0 = x(manfldInd,5);
        vfz0 = x(manfldInd,6);
        
    end
    plot3(x(1:son,1),x(1:son,3),x(1:son,5),'g')
end %this is the section I want to run

% spacecraft initial condition
% x10=0.869581; x30=0; x50=-0.05;
% x20=0.2; x40=-Velo(C,x10,x30,x50); x60=0;
xs0 = -mue;
ys0 = 0.02;
zs0 = 0;
vsx0 = Velo(C,xs0,ys0,zs0)*1.012;
vsy0 = 0;
vsz0 = 0;
%verify zero velocity condition:

% pick a point on the stable manifold to guide to:
% xf0  =  0.4991;
% yf0  = -0.1188;
% zf0  =  0.03149;
% vfx0 =  0;
% vfy0 =  0;
% vfz0 =  0;

%ZEM-ZEV guidance law
% initial conditions
% rs = [xs0;ys0;zs0];
% vs = [vsx0;vsy0;vsz0];
% rf = [xf0;yf0;zf0];
% vf = [vfx0;vfy0;vfz0];
% 
% r  = rf - rs;
% v  = vf - vs;
% 
% simDone = 0;
% dt = 1e-4; % run guidance at 100Hz
% plot3([rs(1) rf(1)],[rs(2) rf(2)],[rs(3) rf(3)],'m--')
% % MAIN LOOP
% tlast = 0;
% tu = .5;
% T = [];
% X = [];
% count = 0;
options = odeset('RelTol',1e-8,'AbsTol',1e-5);
[t,x]=ode45(@cr3bp_3b, [0 4], [xs0 vsx0 ys0 vsy0 zs0 vsz0 ],options);m=length(t);
plot3(x(:,1),x(:,3),x(:,5),'c--');
xGd = x(1,:);
%lookup gravity along trajectory
ti = 0:1e-3:4;
xi(:,1) = interp1(t,x(:,1),ti);
xi(:,2) = interp1(t,x(:,2),ti);
xi(:,3) = interp1(t,x(:,3),ti);
xi(:,4) = interp1(t,x(:,4),ti);
xi(:,5) = interp1(t,x(:,5),ti);
xi(:,6) = interp1(t,x(:,6),ti);

for i = 1:length(ti)
    xp = xi(i,:);
    g(i,:) = gravity(xp);
end


%check to make sure that if I pick something on the stable manifold it
%stays there:
% xf0  =  0.4991;
% yf0  = -0.1188;
% zf0  =  0.03149;
% vfx0 =  Velo(C,xf0,yf0,zf0);
% vfy0 =  0;
% vfz0 =  0;
[t,x]=ode45(@cr3bp_3b, [0 2], [xf0 vfx0 yf0 vfy0 zf0 vfz0 ],options);m=length(t);
plot3(x(:,1),x(:,3),x(:,5),'r*');


% figure(2);
% plot(ti,squeeze(g(:,1)))
cntnue = 1;
dtGd = 5e-2;
time = 0;
tgoGuess = 3;
clear x;
count = 1;
amax = 3;
tol = 1e-4;
Tgo = 2;
while (cntnue)
    
    xs0(count)  = xGd(end,1);
    vsx0(count) = xGd(end,2);
    ys0(count)  = xGd(end,3);
    vsy0(count) = xGd(end,4);
    zs0(count)  = xGd(end,5);
    vsz0(count) = xGd(end,6);
    tme(count)  = time;
    
    sc = [xs0(count) vsx0(count) ys0(count) vsy0(count) zs0(count) vsz0(count) ];
    mp = [xf0 vfx0 yf0 vfy0 zf0 vfz0 ];
    
    norm(sc-mp)
    Tgo = 2-time
    
    if (norm(sc-mp) < tol || Tgo < 0)
        cntnue = 0;
        break;
    end
    [t,x]=ode45(@cr3bp_3b, [0 4], sc,options);m=length(t);
    ti = 0:1e-3:4;
    xi(:,1) = interp1(t,x(:,1),ti);
    xi(:,2) = interp1(t,x(:,2),ti);
    xi(:,3) = interp1(t,x(:,3),ti);
    xi(:,4) = interp1(t,x(:,4),ti);
    xi(:,5) = interp1(t,x(:,5),ti);
    xi(:,6) = interp1(t,x(:,6),ti);
    
    for i = 1:length(ti)
        xp = xi(i,:);
        g(i,:) = gravity(xp);
    end

%     Tgo = fminsearch(@(t) tgo(t,sc,mp,g),tgoGuess)
%     tgoGuess = Tgo;
    gind = min(find(ti>=Tgo));
    G = g(1:gind,:);
    
    acmd = accCmd(sc,mp,G,Tgo);
    norm(acmd);
%     if (norm(acmd)>amax)
%         acmd = acmd/norm(acmd)*amax;
%     end
    
    [t,xGd]=ode45(@fcr3bp_3b, [0 dtGd], sc);m=length(t);
    time = time + dtGd;
    count = count + 1;
    %     for i = 2:length(ti)
    %         f(i) = tgo(ti(i),[xs0 vsx0 ys0 vsy0 zs0 vsz0 ],[xf0 vfx0 yf0 vfy0 zf0 vfz0 ],g);
    %     end
    %     % solve for tgo:
    %     figure(2);
    %     plot(ti,f)
    
     clear x xp g;
end
save trial
plot3(xs0,ys0,zs0,'r')
hold off;
%-----------------------------------
function dx = cr3bp_3b(t, x)
global mue ;
%CR3BP Equations of Motion
%x(1) - x pos
%x(2) - x vel
%x(3) - y pos
%x(4) - y vel
%x(5) - z pos
%x(6) - z vel

dx=zeros(6,1);
r1=sqrt((x(1)+mue)^2+x(3)^2+x(5)^2);
r2=sqrt((x(1)-1+mue)^2+x(3)^2+x(5)^2);
OMx=x(1)-(1-mue)*(x(1)+mue)/r1^3-mue*(x(1)-1+mue)/r2^3;
OMy=x(3)*(1-(1-mue)/r1^3-mue/r2^3);
OMz=-x(5)*((1-mue)/r1^3+mue/r2^3);
dx(1)=x(2);
dx(2)=2*x(4)+OMx;
dx(3)=x(4);
dx(4)=-2*x(2)+OMy;
dx(5)=x(6);
dx(6)=OMz;

function dx = fcr3bp_3b(t, x)
global mue acmd;
%CR3BP Equations of Motion- with commanded acceleration
%x(1) - x pos
%x(2) - x vel
%x(3) - y pos
%x(4) - y vel
%x(5) - z pos
%x(6) - z vel

dx=zeros(6,1);
r1=sqrt((x(1)+mue)^2+x(3)^2+x(5)^2);
r2=sqrt((x(1)-1+mue)^2+x(3)^2+x(5)^2);
OMx=x(1)-(1-mue)*(x(1)+mue)/r1^3-mue*(x(1)-1+mue)/r2^3;
OMy=x(3)*(1-(1-mue)/r1^3-mue/r2^3);
OMz=-x(5)*((1-mue)/r1^3+mue/r2^3);
dx(1)=x(2);
dx(2)=2*x(4)+OMx+acmd(1);
dx(3)=x(4);
dx(4)=-2*x(2)+OMy+acmd(2);
dx(5)=x(6);
dx(6)=OMz+acmd(3);

function f = zvc(x,y)
global mue C
r1=((x + mue).^2+y^2)^(1/2);
r2=((x + mue-1).^2+y^2)^(1/2);
f = 2*((1/2)*(x^2 + y^2)+(1-mue)./r1+mue./r2)-C;

function ind = zerocross(sig)
count = 0;
for i = 1:length(sig)-1
    if ( (sig(i) < 0 && sig(i+1) > 0)  ) %|| ...
         %(sig(i) > 0 && sig(i+1) < 0) )
         count = count + 1;
         ind(count) = i;
    end
end

%---------Jacobi integral--------------
function Vo=Velo(C,xi,yi,zi)
global mue ;
r10=sqrt((xi+mue).^2+yi.^2+zi.^2);r20=sqrt((xi-1+mue).^2+yi.^2+zi.^2);
OM0=0.5*(xi.^2+yi.^2)+(1-mue)/r10+mue/r20+0.5*mue*(1-mue);
Vo=sqrt(2*OM0-C);

%---------tgo function-----------------
function f = tgo(t,x,xf,G)

% if (t<0)
%     t = abs(t);
% end
time = 0:1e-3:4;
vopt = [x(2);x(4);x(6)];
vf = [xf(2);xf(4);xf(6)];
ropt = [x(1);x(3);x(5)];
rf = [xf(1);xf(3);xf(5)];

r = rf - ropt;
v = vf - vopt;

% integrate gravity along trajectory
ind = min(find(time>=t));
g1 = G(1:ind,1);
g2 = G(1:ind,2);
g3 = G(1:ind,3);

Gaccum = zeros(3,1);
try
    Gaccum(1) = trapz(time(1:ind),g1);
    Gaccum(2) = trapz(time(1:ind),g2);
    Gaccum(3) = trapz(time(1:ind),g3);
catch
   fprintf('something\n'); 
end
gt1 = time(1:ind).*g1';
gt2 = time(1:ind).*g2';
gt3 = time(1:ind).*g3';

Gtaccum = zeros(3,1);
Gtaccum(1) = trapz(time(1:ind),gt1);
Gtaccum(2) = trapz(time(1:ind),gt2);
Gtaccum(3) = trapz(time(1:ind),gt3);

pr = 6/t^2*(vopt+vf+Gaccum) - 12/t^3*(r+Gtaccum);
pv = 6/t^2*(r+Gaccum) - 2/t*(2*vf+vopt-Gaccum);

rzem = r - vopt*t-t*Gaccum+Gtaccum;
vzev = v - Gaccum;

g = zeros(3,1);
g = G(1,:)';

f = t^4*(pr'*vopt+pv'*g) ...
  - t^3*(4*pv'*vzev) ...
  + t^2*(4*vzev'*vzev + pv'*rzem) ...
  - t * (12*rzem'*vzev) ...
  +      36*(rzem'*rzem);

d = 1;

%---------lookup gravity---------------
function g = gravity(x)
global mue

g=zeros(3,1);
r1=sqrt((x(1)+mue)^2+x(3)^2+x(5)^2);
r2=sqrt((x(1)-1+mue)^2+x(3)^2+x(5)^2);
OMx=x(1)-(1-mue)*(x(1)+mue)/r1^3-mue*(x(1)-1+mue)/r2^3;
OMy=x(3)*(1-(1-mue)/r1^3-mue/r2^3);
OMz=-x(5)*((1-mue)/r1^3+mue/r2^3);
g(1)=2*x(4)+OMx;
g(2)=-2*x(2)+OMy;
g(3)=OMz;

%compute acceleration command
function acmd = accCmd(x,xf,G,tgo)

time = 0:1e-3:tgo;
vopt = [x(2);x(4);x(6)];
vf = [xf(2);xf(4);xf(6)];
ropt = [x(1);x(3);x(5)];
rf = [xf(1);xf(3);xf(5)];

r = rf - ropt;
v = vf - vopt;

% integrate gravity along trajectory
g1 = G(:,1);
g2 = G(:,2);
g3 = G(:,3);

Gaccum = zeros(3,1);
Gaccum(1) = trapz(time(:),g1(1:length(time)));
Gaccum(2) = trapz(time(:),g2(1:length(time)));
Gaccum(3) = trapz(time(:),g3(1:length(time)));

gt1 = time(:).*g1(1:length(time));
gt2 = time(:).*g2(1:length(time));
gt3 = time(:).*g3(1:length(time));

Gtaccum = zeros(3,1);
Gtaccum(1) = trapz(time(:),gt1);
Gtaccum(2) = trapz(time(:),gt2);
Gtaccum(3) = trapz(time(:),gt3);

rzem = r - vopt*tgo-tgo*Gaccum+Gtaccum;
vzev = v - Gaccum;

acmd = 6/tgo^2*rzem - 2/tgo*vzev;