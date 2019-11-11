function S = vineBeta(d, betaparam)
%// VINE METHOD to generate random correlation matrices
%// with all partial correlations distributed ~ beta(betaparam,betaparam)
%// rescaled to [-1, 1]
%// I found this method on https://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-random-positive-semidefinite-correlation-matrices
%// It is inspired from the paper on https://www.sciencedirect.com/science/article/pii/S0047259X09000876
    P = zeros(d);           %// storing partial correlations
    S = eye(d);

    for k = 1:d-1
        for i = k+1:d
            P(k,i) = betarnd(betaparam,betaparam); %// sampling from beta
            P(k,i) = (P(k,i)-0.5)*2;     %// linearly shifting to [-1, 1]
            p = P(k,i);
            for l = (k-1):-1:1 %// converting partial correlation to raw correlation
                p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
            end
            S(k,i) = p;
            S(i,k) = p;
        end
    end

    %// permuting the variables to make the distribution permutation-invariant
    permutation = randperm(d);
    S = S(permutation, permutation);

return