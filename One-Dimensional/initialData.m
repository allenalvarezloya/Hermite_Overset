function [u1,v1,u2,v2] = initialData(X1,X2,h1,h2,m)

    % Allocate 
    n1 = length(X1);
    n2 = length(X2);
    u1 = zeros(m+1,n1);
    v1 = zeros(m,n1);
    u2 = zeros(m+1,n2);
    v2 = zeros(m,n2);


    % Initial data is u = exp(-20*x^2) & v = 0
    % We map [0,1] to [a,b] by x = (1-r)a + xb 
    % where a = -1.5 and b = 1.5
    % Our transfromation is x = -1.5 + 3r

    for i = 1:n1
        u1(1,i) = exp(-20*(-1.5+3*X1(i))^2);
        u1(2,i) = h1*(180.0 - 360*X1(i))*exp(-20*(-1.5 + 3*X1(i))^2);
        u1(3,i) = (h1^2/factorial(2))*(-360*exp(-20*(-1.5 + 3*X1(i))^2) + (180.0 - 360*X1(i))^2*exp(-20*(-1.5 + 3*X1(i))^2));
    end

    for i = 1:n2
        u2(1,i) = exp(-20*(-1.5+3*X2(i))^2);
        u2(2,i) = h2*(180.0 - 360*X2(i))*exp(-20*(-1.5 + 3*X2(i))^2);
        u2(3,i) = (h2^2/factorial(2))*(-360*exp(-20*(-1.5 + 3*X2(i))^2) + (180.0 - 360*X2(i))^2*exp(-20*(-1.5 + 3*X2(i))^2));
    end
end