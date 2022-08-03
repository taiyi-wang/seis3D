 
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


