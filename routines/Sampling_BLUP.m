function X = Sampling_BLUP(predGP,ybest,xbest,s)
%%%%% construct sampling distribution and sample solution
    d = length(xbest);    %dimension

    %%% Markov Chain Coordinate Sampling (MCCS)
    T = 100;
    X = lhsdesign(s, d);
    for t = 1:T
        Z = X;
        I = (unidrnd(d,s,1)-1)*s + (1:s)'; % Sample Uniformly an interger I from 1 to d
        Z(I) = rand(s,1);  
        [Exz, Vxz] =   predGP([X; Z]);
        Pxz = 1 - normcdf(ybest,Exz,sqrt(Vxz)); % returns the normal cdf at X and Z
        Px = Pxz(1:s);    Pz = Pxz(s+1:end);
        index = rand(s,1)<Pz./Px;
        X(index,:) = Z(index,:);  
    end
    Z = repmat(xbest, s, 1);
    [Exz, Vxz] =   predGP([X; Z]);
    Pxz = 1 - normcdf(ybest,Exz,sqrt(Vxz)); % returns the normal cdf at X and Z
    Px = Pxz(1:s);    Pz = Pxz(s+1:end);
    index = rand(s,1)<Pz./Px;
    X(index,:) = Z(index,:);  
    for t = 1:T
        Z = X;
        I = (unidrnd(d,s,1)-1)*s + (1:s)'; % Sample Uniformly an interger I from 1 to d
        Z(I) = rand(s,1);  
        [Exz, Vxz] =   predGP([X; Z]);
        Pxz = 1 - normcdf(ybest,Exz,sqrt(Vxz)); % returns the normal cdf at X and Z
        Px = Pxz(1:s);    Pz = Pxz(s+1:end);
        index = rand(s,1)<Pz./Px;
        X(index,:) = Z(index,:);  
    end
end
