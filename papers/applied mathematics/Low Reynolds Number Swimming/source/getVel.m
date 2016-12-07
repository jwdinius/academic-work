function velU = getVel(gridStart,...
                       gridEnd  ,...
                       meshWidth,...
                       nr       ,...
                       ns       ,...
                       L0       ,...
                       k        ,...
                       deltar   ,...
                       deltas   ,...
                       xrods    ,...
                       lrod     ,...
                       counter   )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       function getVel                                   %
%       solves Stokes equation at gridpoints defined      %
%       by gridStart,gridEnd, and meshWidth using         %
%       generalized stokeslet and rotlet formulation      %
%                                                         %
%       outputs- velocity field at one timestep           %
%                                                         %
%       Joe Dinius                                        %
%       Adv: Shankar Venkataramani                        %
%       Program in Applied Mathematics                    %
%       25 Nov 2008                                       %
%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
count   = (gridEnd(1) - gridStart(1))/meshWidth + 1;
velU    = zeros(count,count,count,3);
xtr     = zeros(1,3);
pos(1,:)= gridStart(1):meshWidth:gridEnd(1);
pos(2,:)= gridStart(2):meshWidth:gridEnd(2);
pos(3,:)= gridStart(3):meshWidth:gridEnd(3);

%evaluate fluid velocities at grid points
for ii = 1:(floor((gridEnd(1)-gridStart(1))/meshWidth)+1)

    xtr(1) = pos(1,ii);
    for jj = 1:(floor((gridEnd(2)-gridStart(2))/meshWidth)+1)
        xtr(2) = pos(2,jj);
        for kk = 1:(floor((gridEnd(3)-gridStart(3))/meshWidth)+1);
            xtr(3) = pos(3,kk);

            %initialize rotlet and stokeslet value holders
            Ur = [0 0 0];
            Us = [0 0 0];

            %contribution due to rotlet(s) (currently implemented for
            %one motor only)
            for j = 1:nr
                rvec = xtr - xrods(j,:,counter);
                r = norm(rvec);
                Ft= cross(L0(j,:),rvec);
                Ur= Ur + (2*r^2 + 5*deltar^2)*Ft/ ...
                    (16*pi*(r^2 + deltar^2)^(5/2));

            end

            %contribution due to stokeslets
            fjk = 0.0;
            
            for j = ns+1:-1:1
                
                rvec= xtr - xrods(j,:,counter);
                r   = norm(rvec);
                
                if (j == 1)
                    force = fjk; %force at motor
                else
                    xdif= xrods(j,:,counter) - xrods(j-1,:,counter);                    
                    Lkj     = norm(xdif);
                    sprCons = k/(lrod*Lkj) * (Lkj - lrod);
                    fkj     = -[sprCons(1) * xdif(1) ...
                                sprCons(2) * xdif(2) ...
                                sprCons(3) * xdif(3)];
                    force = fkj + fjk; %force at all points connecting two rods
                end
                
                fjk  = -fkj; %utilizing Hooke's law; f_jk = -f_kj

                con1 = (r^2 + 2*deltas^2)/(8*pi*(r^2 + deltas^2)^(3/2));
                con2 =   dot(force,rvec) /(8*pi*(r^2 + deltas^2)^(3/2));

                Us   = Us + (con1*force + con2*rvec);
            end

            velU(ii,jj,kk,:) = Ur + Us;

        end
    end
end