function [tgoEst, cmds, cmdlim] = computeAccCmds(x,tEst)

global finalState controlLimit

rf = finalState([1 3]);
vf = finalState([2 4]);
m  = x(5);
tgoEst = fminunc(@(t0) minimize(t0,x,rf,vf),tEst);

options = odeset('RelTol',1e-8,'AbsTol',1e-5);
[tt,xx]=ode15s(@cr3bp_3b, [0 tgoEst], [x(1) x(2) x(3) x(4) ],options);

ZEM  = rf - xx(end,[1 3]);
ZEV  = vf - xx(end,[2 4]);

acmd = 6/tgoEst^2*ZEM - 2/tgoEst*ZEV;
%acmd = 3/tgo^2*ZEM;

%     ux = acmd(1);
%     uy = acmd(2);
cmdlim = controlLimit / m;
if (norm(acmd) > cmdlim)
    uhat = acmd / norm(acmd);
    ux = cmdlim*uhat(1);
    uy = cmdlim*uhat(2);
else
    ux = acmd(1);
    uy = acmd(2);
end

cmds = [ux;uy];