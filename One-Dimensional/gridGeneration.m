function [X1,X2,h1,h2] = gridGeneration(ratio,m,n,p,q,a,b)
    % For now assume alpha = beta = 0
    
    % The equations x0 + m*h1 = 2*pi & x0 - n*h2 = 0
    % Imply m*h1 + n*h2 = 2*pi => (m+n*ratio)*h1 = 2*pi
    
    h1 = (b-a)/(m+ratio*n);
    h2 = ratio*h1;
    
    % The same equations above imply
    % x0 = (2*pi - (m-ratio*n)*h1)/2
    
    x0 = (b + a - (m - ratio*n)*h1)/2;
    
    X1 = x0 + (0:m)*h1;
    X2 = x0 + (-n:q)*h2;
    
    
    plot(X1,0*X1-1,'x',X2,0*X2+1,'x',X1(1+p),1,'b*',X2(1+n+q),-1,'r*',X1(1),-1,'r*',X2(1+n),1,'b*')
    axis([0 2*pi -2 2])
    

end