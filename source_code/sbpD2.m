function [Dx,A,bcL,bcR] = sbpD2(Nx,dx,cIn,order,x,bcLType,bcRType)
% variable coefficient second derivative operator (compatible, coordinate transform)
% Dx = first derivative operator, d/dx
% A = variable second derivative operator, d(cIn(x) d/dx)/dx, with BC
% bcL = left BC penalty weight
% bcR = right BC penalty weight
%
% Nx = number of grid points
% dx = grid spacing (really, dq when there is a coordinate transform)
% cIn = spatially coefficient matrix
% order = order of accuracy (2 or 4)
% bcLType, bcRtype = boundary condition type (Neumann or Dirichlet)

% 1D matrix factors
xMats = sbp1d_mfc(Nx,dx,order);
xq = diag(xMats.D1int * x);
qx = inv(xq);

% coefficient matrix
c = sparse(diag(cIn(:)));
c = qx.*c;

Dx = qx*xMats.D1int; % 1st derivative

H = xMats.H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create 2nd derivative with variable coefficients
% (This formula comes from equation 8 in
% http://link.springer.com/article/10.1007/s10915-011-9525-z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order==2
  C2x = speye(Nx,Nx);C2x(1,1)=0;C2x(Nx,Nx)=0;
  Rxmu = xMats.D2'*C2x*c*xMats.D2.*(1/4/dx);
elseif order==4
  cc = diag(c);
  B3 = 0.5*sparse(1:Nx,1:Nx,cc,Nx,Nx)...
    + 0.5*sparse(1:Nx,1:Nx,[cc(2:end); cc(end-1)],Nx,Nx);
  
  Rxmu = (1/18/dx).*xMats.D3'*xMats.C3*B3*xMats.D3...
    + (1/144/dx).*xMats.D4'*xMats.C4*c*xMats.D4;
end
Mmux = xMats.D1int'*c*xMats.H*xMats.D1int + Rxmu;
muxBSx = c*xMats.BS;
Dxxmu = xMats.Hinv*(-Mmux + muxBSx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enforcing boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAT penalty weights
if order == 2
  alphaD = -2/dx;
elseif order == 4
  alphaD = -48.0/17.0/dx;
end
alphaT = -1;
beta = 1;

% mapping boundary data to full vector
E0x = sparse(1,1,1,Nx,Nx);
ENx = sparse(Nx,Nx,1,Nx,Nx);
e0x = sparse(1,1,1,Nx,1);
eNx = sparse(Nx,1,1,Nx,1);

if strcmp(bcRType,'Dirichlet'),
  AR = alphaD*c*xMats.Hinv*ENx + beta*xMats.Hinv*(muxBSx)'*ENx;
  bcR = qx*(alphaD*c*xMats.Hinv*eNx + beta*xMats.Hinv*(muxBSx)'*eNx);
elseif strcmp(bcRType,'Neumann'),
  AR = alphaT*xMats.Hinv*ENx*muxBSx;
  bcR = qx*(xMats.Hinv*eNx);
end

if strcmp(bcLType,'Dirichlet'),
  AL = alphaD*c*xMats.Hinv*E0x + beta*xMats.Hinv*(muxBSx)'*E0x;
  bcL = qx*(alphaD*c*xMats.Hinv*e0x + beta*xMats.Hinv*(muxBSx)'*e0x);
elseif strcmp(bcLType,'Neumann'),
  AL = alphaT*xMats.Hinv*E0x*muxBSx;
  bcL = qx*(-xMats.Hinv*e0x);
end

A = qx*(Dxxmu + AL + AR);
end

function mats = sbp1d_mfc(N,h,order)
% creates fully compatible 1D SBP factors for creation of 2D operators


if N == 1,
  if order == 2,
    mats.H = speye(N,N);
    mats.Hinv = sparse(N,N);
    mats.BS = sparse(N,N);
    mats.D1int = sparse(N,N);
    mats.D1 = sparse(N,N);
    mats.D2 = sparse(N,N);
  elseif order == 4
    mats.H = speye(N,N);
    mats.Hinv = sparse(N,N);
    mats.BS = sparse(N,N);
    mats.D1int = sparse(N,N);
    mats.D1 = sparse(N,N);
    mats.D3 = sparse(N,N);
    mats.D4 = sparse(N,N);
    mats.C3 = sparse(N,N);
    mats.C4 = sparse(N,N);
  end
  return
end


switch order
  case(2)
    mats.H = sparse(1:N,1:N,[0.5,ones(1,N-2),0.5]).*h;
    mats.Hinv = sparse(1:N,1:N,[2,ones(1,N-2),2])./h;
    
    % 1st derivative
    mats.D1int = sparse(2:N-1,1:N-2, -0.5*ones(1,N-2),N,N) + ...
      sparse(2:N-1,3:N, 0.5*ones(1,N-2),N,N);
    mats.D1int(1,[1 2]) = [-1 1];
    mats.D1int(end,[end-1 end]) = [-1 1];
    mats.D1int = mats.D1int./h;
    
    mats.BS = sparse(N,N);
    mats.BS(1,:) = -mats.D1int(1,:);
    mats.BS(N,:) = mats.D1int(N,:);
    
    
    % 2nd derivative
    mats.D2 = sparse(2:N-1,2:N-1,-2.*ones(1,N-2),N,N) + ...
      sparse(2:N-1,3:N,ones(1,N-2),N,N) + ...
      sparse(2:N-1,1:N-2,ones(1,N-2),N,N);
    mats.D2(1,1:3) = [1 -2 1];
    mats.D2(N,N-2:N) = [1 -2 1];
    mats.D2 = mats.D2;
    
  case(4)
    vec = [17.0/48.0, 59.0/48.0, 43.0/48.0, 49.0/48.0];
    mats.H = sparse(1:N,1:N,[vec,ones(1,N-8),fliplr(vec)]).*h;
    mats.Hinv = inv(mats.H);
    
    % interior stencil for 1st derivative
    rowV = 5:N-4;
    colV = 3:N-6;
    mats.D1int = sparse(rowV,colV,1/12*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+1,-2/3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+3,2/3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+4,-1/12*ones(size(rowV)),N,N);
    
    mats.D1int(1:4,1:6) = [-24/17,59/34,-4/17,-3/34,0,0; -1/2,0,1/2,0,0,0;...
      4/43,-59/86,0,59/86,-4/43,0; 3/98,0,-59/98,0,32/49,-4/49];
    mats.D1int(N-3:N,N-5:N) = rot90( -mats.D1int(1:4,1:6),2);
    mats.D1int = mats.D1int./h;
    
    % fully compatible
    mats.BS = sparse(N,N);
    mats.BS(1,:) = -mats.D1int(1,:);
    mats.BS(N,:) = mats.D1int(N,:);
    %     mats.BS(1,1:4) = -[-24/17,59/34,-4/17,-3/34];
    %     mats.BS(N,N-3:N) = fliplr([-24/17,59/34,-4/17,-3/34]);
    
    
    % repeating interior stencil [-1, 3, -3, 1]
    rowV = [1,2,4:N-3,N-1,N];
    colV = [1,1,3:N-4,N-3,N-3];
    mats.D3 = sparse(rowV,colV,-1*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+1,3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+2,-3*ones(size(rowV)),N,N)...
      + sparse(rowV,colV+3,1*ones(size(rowV)),N,N);
    % rows 3 and N-3
    vec = [-185893.0/301051.0, 79000249461.0/54642863857.0, -33235054191.0/54642863857.0,...
      -36887526683.0/54642863857.0, 26183621850.0/54642863857.0, -4386.0/181507.0];
    mats.D3(3,1:6) = vec;
    mats.D3(N-2,N-5:N) = -fliplr(vec);
    
    rowV = 1:N;
    colV = [1,1,1,2:N-5,N-4,N-4,N-4];
    mats.D4 = sparse(rowV,colV,1*ones(1,N),N,N)...
      + sparse(rowV,colV+1,-4*ones(1,N),N,N)...
      + sparse(rowV,colV+2,6*ones(1,N),N,N)...
      + sparse(rowV,colV+3,-4*ones(1,N),N,N)...
      + sparse(rowV,colV+4,1*ones(1,N),N,N);
    
    vec1 = [0, 0, 163928591571.0/53268010936.0, 189284.0/185893.0];
    vec2 = [189284.0/185893.0, 0, 163928591571.0/53268010936.0, 0, 0];
    mats.C3 = sparse(1:N,1:N,[vec1,ones(1,N-9),vec2],N,N);
    
    vec = [0, 0, 1644330.0/301051.0, 156114.0/181507.0];
    mats.C4 = sparse(1:N,1:N,[vec,ones(1,N-8),fliplr(vec)],N,N);
    
  otherwise
    display('Error: order can only be 2 or 4.')
    return
end
end