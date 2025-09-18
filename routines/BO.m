function [ output, dmodels] = BO(dmodel, objfun, fitGP, predGP, ACFun, OptACFun, s , maxIter)

    if nargin < 8, maxIter = 30; end
    X = dmodel.X;    Y = dmodel.Y; 
    p = size(X,2);   i = 0;   dmodels = []; 
    while (i < maxIter)
       rng(i)
       [dmodel.ybest, index] = max(Y);    dmodel.xbest = X(index,:);
       if nargout>1, dmodels = [dmodels; dmodel]; end
       fprintf(['Iter: %d; Current Best: ' num2str(dmodel.xbest) ' %f \n'],i, dmodel.ybest);
       switch func2str(ACFun)          
            case 'UCB'
                acHandle = @(x) UCB(x, predGP, dmodel, i+1);
                Xnew = OptACFun(acHandle);
            case 'PI'
                acHandle = @(x) PI(x, predGP, dmodel,dmodel.ybest);
                Xnew = OptACFun(acHandle);
            case 'EI'              
                acHandle = @(x) EI(x, Y, predGP, dmodel);
                Xnew = OptACFun(acHandle);
            case 'Sampling_BLUP'              
                Xnew = Sampling_BLUP(@(X)predGP(dmodel,X), dmodel.ybest, dmodel.xbest, s);
            otherwise
                error('Unknown acquisition function: %s', func2str(ACFun));
       end
       Ynew = objfun(Xnew);
       X = [X; Xnew];   Y = [Y; Ynew];
       t = now; dmodel = fitGP(X, Y, dmodel.theta); dmodel.time = (now-t)*24*3600;
       i = i+1;
    end
    [dmodel.ybest, index] = max(Y);    dmodel.xbest = X(index,:);
    fprintf('Iter: %d; Current Best: [%s] %.3f\n',i, strtrim(num2str(dmodel.xbest, '%.3f ')), dmodel.ybest);
    if nargout>1, dmodels = [dmodels; dmodel]; end
    output = struct('xbest',dmodel.xbest, 'ybest',dmodel.ybest, 'X',X,'Y',Y);
end
