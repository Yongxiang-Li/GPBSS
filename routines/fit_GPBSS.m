function [ fit ] = fit_GPBSS( X, Y, regr, corr, bs, lob, upb, theta0 )
% Periodic Gaussian Process Modeling
% Inputs:
%       X: inputs, which should be normalized in to [0, 1]
%       Y: responses
%       regr: regressor, e.g. @regpoly0
%       corr: correlation function, e.g. @period_sin_gauss_cov
%       lob, upb: lower and upper bounds of the parameter \delta and \theta
%       theta0: initial value for theta without optimizing theta
% Outputs:
%       fit: fitted GPBSS model
    if nargin < 8, theta0 = (lob + upb)/2; end
    if all(X(:) >= 0 & X(:) <= 1)
        X = X;
    else
        X = (X - min(X(:))) / (max(X(:)) - min(X(:)));
    end
  
    for i = 1 : length(bs)
        if bs(i).p<=bs(i).d, bs(i).p = bs(i).d+1; end
    end
    Ut = bs(1).func(X(:,1));
    for i = 2 : size(X,2)
        Ui = bs(i).func(X(:,i));    Ut = sparse_khatrirao(Ut, Ui);
    end
    F = regr(X);
    data = struct('corr',corr, 'regr',regr, 'bs',bs, 'X',X, 'Y',Y, 'F',F, ...
        'Fu',Ut*F, 'Yu',Ut*Y, 'Ut',Ut, 'UtU',full(Ut*Ut'));
    if isempty(lob) && isempty(upb) % without optimization
        [obj, fit] = obj_func(theta0, data);
    else
        [theta, obj, fit, ~] = boxmin(@obj_func, theta0, lob, upb, data);
    end
    
    fit.X = X;    fit.Y = Y;
    fit.theta = theta;
    fit.obj = obj;
    fit.bspline = bs;
    fit.regr = regr;
end

function [obj, fit] = obj_func(para, data) % likelihood at interger
    phi = para(1:end-1); delta = para(end);
    [m, n] = size(data.Ut);
    R1 = data.corr(phi(1), data.bs(1).contrPx);    R1 = R1 + sqrt(eps)*eye(size(R1));
    L = chol(R1, 'lower');    invL = L \ eye(size(R1));    invR = invL'*invL;
    logdetR = 2*m/length(data.bs(1).contrPx)*sum(log(diag(L)));    
    for i = 2 : length(phi)
        Ri = data.corr(phi(i), data.bs(i).contrPx);    Ri = Ri + sqrt(eps)*eye(size(Ri));
        L = chol(Ri, 'lower');    invL = L \ eye(size(Ri));    invR = kron(invR, invL'*invL);
        logdetR = logdetR + 2*m/length(data.bs(i).contrPx)*sum(log(diag(L)));
    end
    [C, flag] = chol(delta^2*invR+data.UtU+sqrt(eps)*eye(m), 'lower');
    if flag~=0
        warning('The covariance matrix is not positive definite')
        obj=inf; return; 
    end
    inv_Fu = C \ data.Fu;    inv_Yu = C \ data.Yu;
    beta = (data.F'*data.F-inv_Fu'*inv_Fu) \ (data.F'*data.Y-inv_Fu'*inv_Yu);
    sigma2 = (norm(data.Y-data.F*beta)^2-norm(inv_Yu-inv_Fu*beta)^2)/(n*delta^2);
    obj = n*log(sigma2) + 2*(n-m)*log(delta) + logdetR + 2*sum(log(diag(C)));
    if ~isreal(obj)
        warning('The knots number is too large')
    end
    if  nargout > 1
        gamma_E = C'\(inv_Yu-inv_Fu*beta);
        R1 = data.corr(phi(1), data.bs(1).contrPx);    R1 = R1 + sqrt(eps)*eye(size(R1));
        R = R1;
        for i = 2 : length(phi)
            Ri = data.corr(phi(i), data.bs(i).contrPx);    Ri = Ri + sqrt(eps)*eye(size(Ri));
            R = kron(R, Ri);
        end
        gamma_V = sigma2*(eye(m) - (C' \ (C \ data.UtU)))*R;
        fit = struct('phi',phi, 'sigma2',sigma2, 'delta',delta, 'beta',beta, ...
            'gamma_E',gamma_E,  'gamma_V',gamma_V);
    end
end

function  [t, f, fit, perf] = boxmin(objfunc, t0, lo, up, data)
%BOXMIN  Minimize with positive box constraints

    % Initialize
    [t, f, itdata] = start(objfunc, t0, lo, up, data);
    if  ~isinf(f)
      % Iterate
      p = length(t);
      if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
      for  k = 1 : kmax
        th = t;
        [t, f, itdata] = explore(objfunc, t, f, itdata, data);
        [t, f, itdata] = move(objfunc, th, t, f, itdata, data);
      end
      [f, fit] = objfunc(t,data);
    end
    perf = struct('nv',itdata.nv, 'perf',itdata.perf(:,1:itdata.nv));
end

function  [t, f, itdata] = start(objfunc, t0, lo, up, data)
% Get starting point and iteration dataameters

    % Initialize
    t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
    D = 2 .^ ([1:p]'/(p+2));
    ee = find(up == lo);  % Equality constraints
    if  ~isempty(ee)
      D(ee) = ones(length(ee),1);   t(ee) = up(ee); 
    end
    ng = find(t < lo | up < t);  % Free starting values
    if  ~isempty(ng)
      t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
    end
    ne = find(D ~= 1);

    % Check starting point and initialize performance info
    f = objfunc(t,data);   nv = 1;
    itdata = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
      'perf',zeros(p+2,200*p), 'nv',1);
    itdata.perf(:,1) = [t; f; 1];
    if  isinf(f)    % Bad dataameter region
      return
    end

    if  length(ng) > 1  % Try to improve starting guess
      d0 = 16;  d1 = 2;   q = length(ng);
      th = t;   fh = f;   jdom = ng(1);  
      for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
          tt = tk .* v; 
          ff = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 1];
          if  ff <= fk 
            tk = tt;  fk = ff;
            if  ff <= f
              t = tt;  f = ff;  jdom = j;
            end
          else
            itdata.perf(end,nv) = -1;   break
          end
        end
      end % improve

      % Update Delta  
      if  jdom > 1
        D([1 jdom]) = D([jdom 1]); 
        itdata.D = D;
      end
    end % free variables

    itdata.nv = nv;
end

function  [t, f, itdata] = explore(objfunc, t, f, itdata, data)
% Explore step

    nv = itdata.nv;   ne = itdata.ne;
    for  k = 1 : length(ne)
      j = ne(k);   tt = t;   DD = itdata.D(j);
      if  t(j) == itdata.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
      elseif  t(j) == itdata.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
      else
        atbd = 0;  tt(j) = min(itdata.up(j), t(j)*DD);
      end
      ff = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 2];
      if  ff < f
        t = tt;  f = ff;
      else
        itdata.perf(end,nv) = -2;
        if  ~atbd  % try decrease
          tt(j) = max(itdata.lo(j), t(j)/DD);
          ff = objfunc(tt,data);  nv = nv+1;
          itdata.perf(:,nv) = [tt; ff; 2];
          if  ff < f
            t = tt;  f = ff;
          else
            itdata.perf(end,nv) = -2;
          end
        end
      end
    end % k

    itdata.nv = nv;
end

function  [t, f, itdata] = move(objfunc, th, t, f, itdata, data)
% Pattern move

    nv = itdata.nv;   ne = itdata.ne;   p = length(t);
    v = t ./ th;
    if  all(v == 1)
      itdata.D = itdata.D([2:p 1]).^.2;
      return
    end

    % Proper move
    rept = 1;
    while  rept
      tt = min(itdata.up, max(itdata.lo, t .* v));  
      ff = objfunc(tt,data);  nv = nv+1;
      itdata.perf(:,nv) = [tt; ff; 3];
      if  ff < f
        t = tt;  f = ff;
        v = v .^ 2;
      else
        itdata.perf(end,nv) = -3;
        rept = 0;
      end
      if  any(tt == itdata.lo | tt == itdata.up), rept = 0; end
    end

    itdata.nv = nv;
    itdata.D = itdata.D([2:p 1]).^.25;
end