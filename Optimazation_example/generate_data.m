function [X] = generate_data(n, p)
    X = lhsdesign(n,p,'criterion','maximin','iterations',1000);
end