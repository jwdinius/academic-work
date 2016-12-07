function vel = interpVel(gridStart,...
                         gridEnd  ,...
                         meshWidth,...
                         ns       ,...
                         velU     ,...
                         xrods    ,...
                         counter   )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       function interpVel                                %
%       interpolate fluid velocity at current location    %
%       of rods                                           %
%                                                         %
%       outputs- flow velocity at rod endpoints           %
%                                                         %
%       Joe Dinius                                        %
%       Adv: Shankar Venkataramani                        %
%       Program in Applied Mathematics                    %
%       25 Nov 2008                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
vel     = zeros(ns+1,3);

for j = 1:ns+1
    
    %interpolate velocity at rod point
    for jj = 1:3
        
        velTemp(:,:,:) = velU(counter,:,:,:,jj);
        vel(j,jj) = interp3(gridStart(1):meshWidth:gridEnd(1),...
                            gridStart(2):meshWidth:gridEnd(2),...
                            gridStart(3):meshWidth:gridEnd(3),...
                            velTemp            ,...
                            xrods(j,1,counter) ,...
                            xrods(j,2,counter) ,...
                            xrods(j,3,counter)  );
                        
    end 

end