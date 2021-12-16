function [u1New,v1New,u2New,v2New] = evolve(h1,h2,dt,m,n,p,q,u1,v1,u2,v2,Hmatu,Hmatv,mx)
% Evolve onto dual grid
u1d = zeros(1+mx,m);
v1d = zeros(mx,m);
u2d = zeros(1+mx,n+q);
v2d = zeros(mx,n+q);
for i = 1:m
    interp_u1 = Hmatu*[u1(:,i); u1(:,i+1)];
    interp_v1 = Hmatv*[v1(:,i); v1(:,i+1)];
    [u,v] = updateCoeffs(interp_u1,interp_v1,mx,h1,dt);
    u1d(:,i) = u(1:mx+1);
    v1d(:,i) = v(1:mx);
end
for i = 1:n+q
    interp_u2 = Hmatu*[u2(:,i); u2(:,i+1)];
    interp_v2 = Hmatv*[v2(:,i); v2(:,i+1)];
    [u,v] = updateCoeffs(interp_u2,interp_v2,mx,h2,dt);
    u2d(:,i) = u(1:mx+1);
    v2d(:,i) = v(1:mx);
end
% Evolve back to primal
for i = 1:m-1
    interp_u1 = Hmatu*[u1d(:,i); u1d(:,i+1)];
    interp_v1 = Hmatv*[v1d(:,i); v1d(:,i+1)];
    [u,v] = updateCoeffs(interp_u1,interp_v1,mx,h1,dt);
    u1(:,1+i) = u(1:mx+1);
    v1(:,1+i) = v(1:mx);
end
for i = 1:n+q-1
    interp_u2 = Hmatu*[u2d(:,i); u2d(:,i+1)];
    interp_v2 = Hmatv*[v2d(:,i); v2d(:,i+1)];
    [u,v] = updateCoeffs(interp_u2,interp_v2,mx,h2,dt);
    u2(:,1+i) = u(1:mx+1);
    v2(:,1+i) = v(1:mx);
end

% Inject data from u1 to u2 and u2 to u1
for idx = 0:mx
    u2(1+idx,1+n+q) = (h2/h1)^idx*u1(1+idx,1+p);
    u1(1+idx,1) = (h1/h2)^idx*u2(1+idx,1+n);
end
for idx = 0:mx-1
    v2(1+idx,1+n+q) = (h2/h1)^idx*v1(1+idx,1+p);
    v1(1+idx,1) = (h1/h2)^idx*v2(1+idx,1+n);
end
% Reflection Boundary conditions
ughost = zeros(mx+1,1);
vghost = zeros(mx,1);
% Right
for idx = 0:mx
    ughost(1+idx) = (-1)^(1+idx)*u1d(1+idx,m);
end
for idx = 0:mx-1
    vghost(1+idx) = (-1)^(1+idx)*v1d(1+idx,m);
end
interp_u = Hmatu*[u1d(:,m); ughost];
interp_v = Hmatv*[v1d(:,m); vghost];

[u,v] = updateCoeffs(interp_u,interp_v,mx,h1,dt);
u1(:,1+m) = u(1:mx+1);
v1(:,1+m) = v(1:mx);
% Left
for idx = 0:mx
    ughost(1+idx) = (-1)^(1+idx)*u2d(1+idx,1);
end
for idx = 0:mx-1
    vghost(1+idx) = (-1)^(1+idx)*v2d(1+idx,1);
end
interp_u = Hmatu*[ughost; u2d(:,1)];
interp_v = Hmatv*[vghost; v2d(:,1)];

[u,v] = updateCoeffs(interp_u,interp_v,mx,h2,dt);
u2(:,1) = u(1:mx+1);
v2(:,1) = v(1:mx);


    
u1New = u1;
v1New = v1;
u2New = u2;
v2New = v2;
end