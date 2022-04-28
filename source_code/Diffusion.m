classdef Diffusion

% solves linear diffusion equation of form:
% B*dp/dt + dq/dx = 0, q = -K*dp/dx,
% B = b*A, where b=compressibility times porosity, A=cross-sectional area
% K = k*A/mu, where k=permeability, mu=viscosity, A=cross-sectional area
% Note that B and K can depend on x.
% Also, x and p are stored externally and passed to methods in Diffusion class.
    
    properties
        order_diffusion = 4; % FD order of accuracy
        nx % number of grid points
        B % coefficient in diffusion PDE, stored as sparse diagonal nx by nx matrix
        K % coefficient in diffusion PDE, stored as sparse diagonal nx by nx matrix
        Sw % well storage capacity (connected at left boundary)
        BC % right/outer boundary condition (pressure or flowrate)
           % below are used to solve PDE, not input by user
        D1 % first derivative: d/dx
        D2 % variable coefficient second derivative: d/dx (B d/dx)
        s % SAT penalty vector for enforcing left boundary condition
        S % matrix version of s
    end
    
    
    methods
        
        
        function D = Diffusion(nx,B,K,Sw,BC)
        % Diffusion constructor
            D.nx = nx;
            D.B = B;
            D.K = K;
            D.Sw = Sw;
            D.BC = BC;
        end
        
        
        function D = setParam(D,param,value)
        % set parameter for diffusion object
            switch param
              case 'nx'
                D.nx = value;
              case 'B'
                D.B = value;
              case 'K'
                D.K = value;
              case 'Sw'
                D.Sw = value;
              case 'BC'
                D.BC = value;
              otherwise
                error('invalid diffusion object parameter to set')
            end
        end
        
        
        function value = getParam(P,param)
	% get parameter for diffusion object
            switch param
              case 'nx'
                value = D.nx;
              case 'B'
                value = D.B;
              case 'K'
                value = D.K;
              case 'Sw'
                value = D.Sw;
              case 'BC'
                value = D.BC;
              otherwise
                error('invalid dffusion object parameter')
            end
        end
        
        
        function D = discretize(D,x)
        % finite difference discretization of diffusion equation
            
            assert(length(x)==D.nx,'length(x) must match nx in Diffusion class')
            
            % SBP-SAT operators
            
            % B*dp/dt + dq/dx = 0, q = -K*dp/dx
            % convert B,K to sparse diagonal matrices if they are not already
            
            [n1,n2] = size(D.B);
            if n1==D.nx&n2==1
                D.B = spdiags(D.B,0,D.nx,D.nx);
            elseif n1==1&n2==D.nx
                D.B = spdiags(D.B',0,D.nx,D.nx);
            elseif n1==D.nx&n2==D.nx
                assert(isdiag(D.B),'B must be diagonal in Diffusion class')
                D.B=sparse(D.B);
            else
                error('B has invalid size')
            end
            
            [n1,n2] = size(D.K);
            if n1==D.nx&n2==1
                D.K = spdiags(D.K,0,D.nx,D.nx);
            elseif n1==1&n2==D.nx
                D.K = spdiags(D.K',0,D.nx,D.nx);
            elseif n1==D.nx&n2==D.nx
                assert(isdiag(D.K),'K must be diagonal in Diffusion class')
                D.K=sparse(D.K);
            else
                error('K has invalid size')
            end
            
            D = D.calcD1D2(diag(D.K),x); % calculate first and second derivative operators (see function below)
            
        end
        
        
        function D = calcD1D2(D,K,x)
        % calculate first derivative D1 and variable coefficient second derivative operator D2
        % K = variable coefficient in d/dx (K d/dx) operator
        % x = distance
            
            h = (max(x)-min(x))/(D.nx-1); % grid spacing (for uniform grid, otherwise for x=x(q), h=dq)
            
            switch D.BC
              case 'pressure'
                [D1,D2,s,~]=sbpD2(D.nx,h,diag(D.K),D.order_diffusion,x,'Neumann','Dirichlet');
              case 'flowrate'
                [D1,D2,s,~]=sbpD2(D.nx,h,diag(D.K),D.order_diffusion,x,'Neumann','Neumann');
              otherwise
                error('invalid BC in Diffusion object')
            end
            D.D1 = sparse(D1); D.D2 = sparse(D2);
            D.s = -s; D.S = sparse(D.nx,D.nx); D.S(:,1) = D.s;
            
        end
          
        
        function p = update_diffusion(D,p,dt,Q)
        % backward Euler update to diffusion equation
            
        % INPUT:
        % D = diffusion class object
        % p = pressure at previous time step at all points in object
        % dt = time step
        % Q = flow rate in at left boundary
            
        % OUTPUT:
        % p = pressure at current time step, after update (this gets overwritten from input)

        p = (D.B+D.S*D.Sw-dt*D.D2)\(D.B*p+dt*D.s*Q);
            
        end
        
        function q = flow_rate(D,p)
        % calculate flow rate at left/inner boundary
            q = -D.K(1,1)*D.D1(1,:)*p;
        end   
        
        
        function q = flow_rate_everywhere(D,p)
        % calculate flow rate
            q = -D.K*D.D1*p;
        end   
        
     
    end
end