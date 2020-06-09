clc
clear
clear all
close all
clc

k = 1e-2;    % time step
tf = 20;    % final time
nx = 30;      % number of x-gridpoints
ny = 30;      % number of y-gridpoints

h=1/(nx-1); %spatial step
nt=0.1;   %time step for plot
dnt=nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 (ALREADY IMPLEMENTED. DO NOT MODIFY)
%
%
x=linspace(0,1,nx);
y=linspace(0,1,ny);
[X,Y]=meshgrid(x,y);
X=X(end:-1:1,:);
Y=Y(end:-1:1,:);
u=exp(-100*((X-0.5).^2+(Y-0.5).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% STEP 2
for t=0:k:tf
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% STEP 2a
%
for i=2:ny-1
    for j=2:nx-1
        
u(i,j)=u(i,j)+k*(-0.02*(u(i,j+1)-u(i,j-1)+u(i-1,j)-u(i+1,j))/(2*h))+k*(1/2000)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-4*u(i,j))/(h^2); 

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% plot section (DO NOT MODIFY)
    if t>=nt-1e-4
        clf
        [Xi,Yi]=meshgrid(linspace(0,1,500),linspace(0,1,500));
        ui=interp2(X,Y,u,Xi,Yi);
        surf(Xi,Yi,ui,'EdgeColor','none');
        zlim([0,1])
        hold on
        title(['time: t=',num2str(t)],'FontSize',20)
        pause(1e-5)
        nt=nt+dnt;
    end
end
