function F = minimize(tt,x,rf,vf)
options = odeset('RelTol',1e-8,'AbsTol',1e-5);
[ttt,xxx]=ode15s(@cr3bp_3b, [0 tt], [x(1) x(2) x(3) x(4)],options);

vec1 = [rf(1)-xxx(end,1);...
    vf(1)-xxx(end,2);...
    rf(2)-xxx(end,3);...
    vf(2)-xxx(end,4)];

F = vec1'*vec1;
end