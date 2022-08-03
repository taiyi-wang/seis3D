function [ta, Va, Da, Psia, Pa, taua] = springslider(M, D0, Psi0, tauy0, tmax, t_s, p_s, dtaux_s, dtauy_s, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive time stepping to solve spring-slider ODE system

% Input:
% M = input structure M1 with f0, V0, a, b, d_c, eta, tauy_0
% D0 = initial displacement
% Psi0 = initial state variable
% tauy0 = initial normal stress 
% tmax = maximum time for calculation
% t_s = Nstep X 1 vector of time (s) from the main fault only calculation
% p_s = Nstep X 1 vector of pore pressure (Pa) from the main fault only calculation
% dtaux_s = Nstep X 1 vector of shear stress change in x (along dip, Pa) 
% dtauy_s = Nstep X 1 vector of normal stress change in y (fault normal, Pa)
% k       = scalar stiffness relating slip on spring slider with shear stress change on spring slider (Pa/m)

% Output:
% ta = Mstep X 1 vector of time (s) from the adaptive time stepping
% Va = Mstep X 1 vector of velocities (m/s) 
% Da = Mstep X 1 vector of displacement (m)
% Psia = Nstep X 1 vector of state variable
% Pa = Mstep x 1 vector of pore pressure at spring slider (Pa)
% taua = Nstep X 1 vector of shear stress on spring slider (Pa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% using adaptive time-stepping
tol = M.tol; % tolerance
dt = M.dt; % initial time step (s)
dtmax = M.dtmax; % maximum time step (s)
safety = M.safety; % safety factor
   
t=0;
    
% store solution
ta=t; Da=D0; Psia=Psi0;

tauya = tauy0 + interp1(t_s, dtauy_s, t);
pa = interp1(t_s, p_s, t);
dtauxa = interp1(t_s, dtaux_s, t);

[V1,G1,tau1] = sliderODE(Da(end),Psia(end),t, dtauxa, tauya-pa, k, M);
Va=V1; taua=tau1; % stage 1 values are stored
Pa=pa;
    
err=0; dta=dt;

tStart = tic; % record computation time
while t<tmax && toc(tStart)<60

    % adjust dt to stop at tmax
    if t+dt>tmax, dt=tmax-t; end
        
    % three-stage method with embedded error estimate
    % Note: tayta, pa, and dtauxa have to be updated at intermediate time
    % steps because they feedback into evaluation of V and G at 
    % intermediate time steps
    tauya = tauy0 + interp1(t_s, dtauy_s, t+0.5*dt);
    pa = interp1(t_s, p_s, t+0.5*dt);
    dtauxa = interp1(t_s, dtaux_s, t+0.5*dt);
    [V2,G2] = sliderODE(Da(end)+0.5*dt*V1,Psia(end)+0.5*dt*G1,t+0.5*dt, dtauxa, tauya-pa, k, M);
    
    tauya = tauy0 + interp1(t_s, dtauy_s, t+dt);
    pa = interp1(t_s, p_s, t+dt);
    dtauxa = interp1(t_s, dtaux_s, t+dt);
    [V3,G3] = sliderODE(Da(end)+dt*(-V1+2*V2),Psia(end)+dt*(-G1+2*G2),t+dt, dtauxa, tauya-pa, k,  M);
        
    % second order update
    D2   = Da(end)  +dt/2*(V1+V3);
    Psi2 = Psia(end)+dt/2*(G1+G3);
        
    % third order update
    D3   = Da(end)  +dt/6*(V1+4*V2+V3);
    Psi3 = Psia(end)+dt/6*(G1+4*G2+G3);
        
    q = 2; % order of accuracy of lower order update
        
    % local error estimate
    er = norm([D2-D3; Psi2-Psi3]);
        
    if er<tol
        % update solution
        t = t+dt; ta=[ta; t];
        Da = [Da; D3]; Psia = [Psia; Psi3]; % use third-order update
            
        % store error and time step
        err=[err; er]; dta=[dta; dt];
            
        % evaluate stage 1 values for next time step
        tauya = tauy0 + interp1(t_s, dtauy_s, t);
        pa = interp1(t_s, p_s, t);
        dtauxa = interp1(t_s, dtaux_s, t);
    
        [V1,G1,tau1] = sliderODE(Da(end),Psia(end),t, dtauxa, tauya-pa, k, M);
        Va=[Va; V1]; taua=[taua; tau1]; % stage 1 values are stored
        Pa=[Pa; pa];
    end
        
    % adjust time step
    dt = safety*dt*(tol/er)^(1/(q+1));
    dt = min(dt,dtmax);
        
end
end

% functions below here

function [V,G,tau] = sliderODE(D,Psi,t, dtaux, tauy, k, M)
    
    % evaluate stress when V=0

    tauLock = M.tau0+dtaux+k*D;

    % set bounds on V for root-finding

    if tauLock>0
        Vmin = 0; Vmax = tauLock/M.eta;
    else
        Vmin = tauLock/M.eta; Vmax = 0;
    end

    % solve stress=strength for V
    
    atol = 1e-14; rtol = 1e-6;
    
    V = hybrid(@(V) solveV_ss(V,tauLock, tauy, Psi,M) ,Vmin,Vmax,atol,rtol);
    
    % then evaluate total shear stress

    tau = tauLock-M.eta*V;

    % and state evolution, G = dPsi/dt

    if V==0
        G = 0; % special case to avoid log(0)
    else
        % use aging law
        G = (M.b.*M.V0./M.dc).*(exp((M.f0 - Psi)./M.b) - V./M.V0);
    end
        
end


function residual = solveV_ss(V,tauLock, tauy, Psi,M)

    stress = tauLock-M.eta*V;
    f = M.a*asinh(V/(2*M.V0)*exp(Psi/M.a));
    strength = f*tauy;
    residual = stress-strength;
    
end


function [x,err]=hybrid(func,a,b,atol,rtol)

  % hybrid method solves func(x)=0 for some root x within (a,b)
  % returns x, estimate of root with absolute error less than atol
  % or relative error less than rtol
  
  % function values at endpoints
  fa = func(a);
  fb = func(b);

  % make sure root is bracketed; otherwise return
  if sign(fa)==sign(fb) | isnan(fa) | isnan(fb)
    disp('error: root not bracketed or function is NaN at endpoint')
    x = NaN; err = NaN;
    return
  end

  % set up secant method, storing old values as xold and fold, new ones as x and f
  % use bisection brackets to start secant (this is somewhat arbitrary)
  xold = a;
  fold = fa;
  x = b;
  f = fb;
  
  % begin iterations,
  % keeping track of error at each iteration in vector err
  n = 0; err = [];
  update = 'input'; % character string stating type of update used in previous interation
  while b-a>atol+rtol*abs(x) % safe to have infinite loop since bisection guaranteed to converge
      
      err = [err b-a]; % add to end of vector the current error (interval width)

      % formatted printing so you can watch method converge
      %fprintf('%6i %20.10f %20.10f %20.10f %s\n',n,a,x,b,update)

      n = n+1; % iteration number

      % first calculate (tenative) secant update
      dfdx = (f-fold)/(x-xold); % approximation to df/dx
      dx = -f/dfdx; % update to x
      xs = x+dx; % secant update
      
      % determine if secant method will be used
      if (xs<a) | (xs>b)
          use_secant = false;  % not if update outside (a,b)
      else
          fs = func(xs); % function value at secant update
          % calculate interval reduction factor = (old interval width)/(new interval width)
          if sign(fs)==sign(fa)
              IRF = (b-a)/(b-xs); % would update a=xs
          else
              IRF = (b-a)/(xs-a); % would update b=xs
          end
          if IRF<2
              use_secant = false;
          else
              use_secant = true;
          end
      end

      xold = x; fold = f; % store these values for next iteration

      % now update
      if use_secant
          update = 'secant';
          x = xs;
          f = fs;
      else
          update = 'bisection';
          x = (a+b)/2; % midpoint
          f = func(x); % function value at midpoint
      end
      
      % update one endpoint based on sign of function value at updated x
      if sign(f)==sign(fa)
          a = x;
      else
          b = x;
      end
      
  end
  
end

