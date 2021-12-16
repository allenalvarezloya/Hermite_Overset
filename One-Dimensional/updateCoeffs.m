function [unew, vnew] = updateCoeffs(u,v,m,h,dt)
U = zeros(2*m+5,2*m+2);
V = zeros(2*m+4,2*m);
U(1:2*m+2,1) = u;
V(1:2*m,1) = v;

% Recurstion relation
for s = 1:2*m+1
    for l = 0:2*m+2
        U(1+l,1+s) = (dt/s)*V(1+l,1+s-1);
        V(1+l,1+s) = ((l+2)*(l+1)/s)*(dt/h^2)*U(1+l+2,1+s-1);
    end
end
% Timestep
timeStep = (0.5).^(0:2*m+1);
unew = zeros(m+1,1);
vnew = zeros(m,1);
% Update u
for s = 0:2*m+1
    unew(1+s,1) = sum(timeStep.*U(1+s,:));
end
for s = 0:2*m-1
    vnew(1+s,1) = sum(timeStep.*V(1+s,:));
end


end