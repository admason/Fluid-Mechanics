%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Re = 1000;     % Reynolds number (CHANGE ACCORDINGLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1e-3;    % time step
tf = 15;    % final time
nx = 40;      % number of x-gridpoints
ny = 40;      % number of y-gridpoints
 
 
h=1/(nx-1); %spatial step
uTOP=1;   %velocity of the lid
nt=0.1;   %time step for plot
dnt=nt;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 (ALREADY IMPLEMENTED. DO NOT MODIFY)
%
%
div_U_star=zeros(ny,nx);
u=zeros(ny,nx);
u(1,:)=uTOP;
v=zeros(ny,nx);
ustar=zeros(ny,nx);
vstar=zeros(ny,nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
[~,M_P]=PoissonSolver(div_U_star,h);
for t=0:k:tf
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2a
%
%
for i=2:ny-1
    for j=2:nx-1
        ustar(i,j)=u(i,j)+k*(-((((u(i,j+1))^2)-((u(i,j-1))^2)+((u(i-1,j))*(v(i-1,j)))-((u(i+1,j))*(v(i+1,j))))/(2*h))+((1/Re)*((u(i-1,j))+(u(i+1,j))+(u(i,j-1))+(u(i,j+1))-(4*(u(i,j))))/(h^2)));
        vstar(i,j)=v(i,j)+k*(-((((u(i,j+1))*(v(i,j+1)))-(((u(i,j-1))*(v(i,j-1))))+((v(i-1,j))^2)-((v(i+1,j))^2))/(2*h))+((1/Re)*((v(i-1,j))+(v(i+1,j))+(v(i,j-1))+(v(i,j+1))-(4*(v(i,j))))/(h^2)));
    end
end
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2b
%
%
 
for i=2:ny-1
    for j=2:nx-1
        div_U_star(i,j)=(((ustar(i,j+1))-(ustar(i,j-1))+(vstar(i-1,j))-(vstar(i+1,j)))/(2*h));
    end
end
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2c (ALREADY IMPLEMENTED. DO NOT MODIFY)
    pp = M_P\[div_U_star(:);0]; p = reshape(pp(1:end-1),ny,nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2d 
%    
%
for i=2:ny-1
    for j=2:nx-1
        u(i,j)=ustar(i,j)-(((p(i,j+1))-(p(i,j-1)))/(2*h));
        v(i,j)=vstar(i,j)-(((p(i-1,j))-(p(i+1,j)))/(2*h));
    end
end
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot section (DO NOT MODIFY)
    if t>=nt-1e-4
        clf
        x=linspace(0,1,nx);
        y=linspace(0,1,ny);
        [X,Y]=meshgrid(x,y);
        X=X(end:-1:1,:);
        Y=Y(end:-1:1,:);
        [Xi,Yi]=meshgrid(linspace(0,1,300),linspace(0,1,300));
        magnV=sqrt(u.^2+v.^2);
        magnVi=interp2(X,Y,magnV,Xi,Yi);
        surf(Xi,Yi,magnVi-max(magnVi(:))-1,'EdgeColor','none');
        hold on
        step=1;
        quiver(X(1:step:end,1:step:end),Y(1:step:end,1:step:end),u(1:step:end,1:step:end)./sqrt(u(1:step:end,1:step:end).^2+v(1:step:end,1:step:end).^2+eps),v(1:step:end,1:step:end)./sqrt(u(1:step:end,1:step:end).^2+v(1:step:end,1:step:end).^2+eps),'w-')
        hold off, axis equal, axis([0 1 0 1])
        title(['time: t=',num2str(t)],'FontSize',20)
        pause(1e-5)
        nt=nt+dnt;
    end
end
