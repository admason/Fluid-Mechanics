function [p,M,RHS] = PoissonSolver(RHS_Lap_p,h)

% global BC driven_cavity;
dx=h;
dy=h;
Nx=size(RHS_Lap_p,2)-1;
Ny=size(RHS_Lap_p,1)-1;
nn = Ny+1;
mm = Nx+1;
% M = speye(nn*mm+1);
RHS = zeros(nn*mm+1,1);
nnz = (5+2)*(Nx-1)*(Ny-1)+8*(Nx-1)+8*(Ny-1)+4;
rows = zeros(nnz,1);
cols = zeros(nnz,1);
vals = zeros(nnz,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont = 0;
for i = 2:Ny
    for j = 2:Nx
        k = i+(j-1)*nn;
        % STO SUPPONENDO dx=dy !!!
        %M(k,[k k-nn k+nn k-1 k+1])=[4 -1 -1 -1 -1]/dx^2;
        rows(cont+(1:5)) = [k k k k k];
        cols(cont+(1:5)) = [k k-nn k+nn k-1 k+1];
        vals(cont+(1:5)) = [4 -1 -1 -1 -1]/dx^2;
        cont = cont + 5;
        RHS(k)=RHS_Lap_p(k);
%         M(nn*mm+1,k)=1;
%         M(k,nn*mm+1)=dx;
        rows(cont+(1:2)) = [nn*mm+1 k];
        cols(cont+(1:2)) = [k nn*mm+1];
        vals(cont+(1:2)) = [1 dx];
        cont = cont + 2;
        
        
    end
end


for j = 2:Nx
    i = 1;
    k = i+(j-1)*nn;
    %M(k,[k k+1 k+2])=[3 -4 1]/2/dy;
    rows(cont+(1:3)) = [k k k];
    cols(cont+(1:3)) = [k k+1 k+2];
    vals(cont+(1:3)) = [3 -4 1]/2/dy;
    cont = cont + 3;
        
    RHS(k) = 0;
    %M(k,nn*mm+1)=1;
    rows(cont+1) = k;
    cols(cont+1) = nn*mm+1;
    vals(cont+1) = 1;
    cont = cont + 1;
    
    i=nn;
    k = i+(j-1)*nn;
%     M(k,[k k-1 k-2])=[3 -4 1]/2/dy;
    rows(cont+(1:3)) = [k k k];
    cols(cont+(1:3)) = [k k-1 k-2];
    vals(cont+(1:3)) = [3 -4 1]/2/dy;
    cont = cont + 3;
    
    RHS(k)=0;
%     M(k,nn*mm+1)=1;
    rows(cont+1) = k;
    cols(cont+1) = nn*mm+1;
    vals(cont+1) = 1;
    cont = cont + 1;
end
    
for i = 2:Ny
    j = 1;
    k = i+(j-1)*nn;
%     M(k,[k k+nn k+2*nn]) = [3 -4 1]/2/dx;
    rows(cont+(1:3)) = [k k k];
    cols(cont+(1:3)) = [k k+nn k+2*nn];
    vals(cont+(1:3)) = [3 -4 1]/2/dx;
    cont = cont + 3;
    
    RHS(k) = 0;
%     M(k,nn*mm+1) = 1;
    rows(cont+1) = k;
    cols(cont+1) = nn*mm+1;
    vals(cont+1) = 1;
    cont = cont + 1;

    j = mm;
    k = i+(j-1)*nn;
%     M(k,[k k-nn k-2*nn]) = [3 -4 1]/2/dx;
    rows(cont+(1:3)) = [k k k];
    cols(cont+(1:3)) = [k k-nn k-2*nn];
    vals(cont+(1:3)) = [3 -4 1]/2/dx;
    cont = cont + 3;
    
    RHS(k) = 0;
%     M(k,nn*mm+1) = 1;
    rows(cont+1) = k;
    cols(cont+1) = nn*mm+1;
    vals(cont+1) = 1;
    cont = cont + 1;

end

RHS(nn*mm+1) = 0;
rows(cont+(1:4)) = [1 nn nn*mm-nn+1 nn*mm];
cols(cont+(1:4)) = [1 nn nn*mm-nn+1 nn*mm];
vals(cont+(1:4)) = 1; %oppure [1 1 1 1] lo capisce da solo

M = -sparse(rows,cols,vals);
press = M\RHS;

p = reshape(press(1:nn*mm),nn,mm);
