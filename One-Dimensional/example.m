z = 3;  % Number of refinements
a = 0.0; 
b = 1.0;
m = 32*z;
n = 8*z;
ratio = 0.5;
p = 5;
q = 10;
mx = 3;
Hmatu = Hermite_map(mx,0,1,0.5,0);
Hmatv = Hermite_map(mx-1,0,1,0.5,0);
[X1,X2,h1,h2] = gridGeneration(ratio,m,n,p,q,a,b);
[u1,v1,u2,v2] = initialData(X1,X2,h1,h2,mx);
dt = 0.9*min(h1,h2);
T = 2.0;
t = 0.0;
while t < T
    if(t + dt > T)
        dt = T - t;
    end
    [u1,v1,u2,v2] = evolve(h1,h2,dt,m,n,p,q,u1,v1,u2,v2,Hmatu,Hmatv,mx);
    t = t + dt;
    plot(X1,u1(1,:),'b',X2,u2(1,:),'r')
    axis([0 1 -1 1])
    TITLE = sprintf('t = %0.3f',t);
    title(TITLE)
    drawnow
end





