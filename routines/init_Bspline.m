function [ bspline ] = init_Bspline(d, Nk)
    bspline = [];
    for i = 1 : length(Nk)
        bs.d = d;    bs.p = Nk(i);
        bs.knots = [zeros(1,bs.d) equalspace(0,1,bs.p-bs.d+1) ones(1,bs.d)]';
        bs.contrPx = equalspace(0,1,bs.p)';
        bs.func = @(x)BsplineMatrix(bs.d, bs.p, bs.knots, x);
        bspline = [bspline; bs];
    end
end