function [Ye, Yv ] = predict_GPBSS(model, X)
    F = model.regr(X); 
    Ut = model.bspline(1).func(X(:,1));
    for i = 2 : size(X,2)
        Ui = model.bspline(i).func(X(:,i));    Ut = sparse_khatrirao(Ut, Ui);
    end   
    Ye = F*model.beta+Ut'*model.gamma_E;
    if nargout==2
        Yv = full(dot(Ut, (model.gamma_V*Ut)))';
    end
end