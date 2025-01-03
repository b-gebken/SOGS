function y = halfhalf(x)

    A = diag([1,0,1,0,1,0,1,0]);
    B = diag(1./(1:8).^2);

    y = sqrt(x'*A*x) + x'*B*x;

end

