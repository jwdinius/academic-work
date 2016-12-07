%plotResults
close all; clear all; clc

load dumbbell;
count = 1;

[x,y,z] = ndgrid(0:.1:1,0:.1:1,0:.1:1);
u = zeros(11,11,11);
v = zeros(11,11,11);
w = zeros(11,11,11);
%velU = 100*velU;
for i = 1:length(velU(:,1,1,1,1));
    
    for ii=1:11
        for jj = 1:11
            for kk = 1:11
                u(ii,jj,kk) = velU(i,ii,jj,kk,1);
                v(ii,jj,kk) = velU(i,ii,jj,kk,2);
                w(ii,jj,kk) = velU(i,ii,jj,kk,3);
            end
        end
    end
    
    figure(1)
    
%     quiver3(x,y,z,u,v,w);
%     xlim([0 1]);
%     ylim([0 1]);
%     zlim([0 1]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     title(['t = ' num2str(i*.01) ' sec']);
    
    for j = 3:3
        figure(1)
        xslice(:,:) = x(:,j,:);
        yslice(:,:) = y(j,:,:);
        zslice(:,:) = z(:,j,:);
        z1slice(:,:) = z(j,:,:);
        uslice(:,:) = u(:,j,:);
        vslice(:,:) = v(j,:,:);
        wslice(:,:) = w(:,j,:);
        w1slice(:,:)= w(j,:,:);
%         quiver(xslice,zslice,uslice,wslice);
%         axis([-.2 1.2 -.2 1.2])
%         xlabel('x')
%         ylabel('z')
%        title(['y = ' num2str(.1*(j-1))]);
%        legend(['t = ' num2str(i)]);
%         figure(2);
        quiver(yslice,z1slice,vslice,w1slice);
        xlabel('y')
        ylabel('z')
        title(['x = ' num2str(.1*(j-1))]);
        axis([-.2 1.2 -.2 1.2])        
        legend(['t = ' num2str(i)]);
   end
    M(i) = getframe;
%     figure(100);
%     plot3(xrods(:,1,i),xrods(:,2,i),xrods(:,3,i));
%     xlim([-1 1]);
%     ylim([-1 1]);
%     zlim([-1 1]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     title(['t = ' num2str(i*.01) ' sec']);
    
    count = count + 1;
    
end
movie(M);
movie2avi(M,'dumbbell.avi');