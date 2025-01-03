function [y,I] = mxhilb(x)

    n = size(x,1);
    tmp = zeros(n,1);
    for i = 1:n
        tmp(i) = sum(x./(i + (1:n)' - 1),1);
    end

    [y,I] = max(abs(tmp),[],1);
    s = sign(tmp(I(1)));
    if(s == 0)
        s = 1;
    end
    I = s*I;
    
end

