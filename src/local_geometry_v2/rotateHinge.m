function M = rotateHinge(u,phi)

u = u/norm(u);
M = zeros(3);

M(1,1) = cos(phi) + (1-cos(phi))*u(1)^2;
M(1,2) = (1-cos(phi))*u(1)*u(2) - u(3)*sin(phi);
M(1,3) = (1-cos(phi))*u(1)*u(3) + u(2)*sin(phi);

M(2,1) = (1-cos(phi))*u(1)*u(2) + u(3)*sin(phi);
M(2,2) = cos(phi) + (1-cos(phi))*u(2)^2;
M(2,3) = (1-cos(phi))*u(2)*u(3) - u(1)*sin(phi);

M(3,1) = (1-cos(phi))*u(1)*u(3) - u(2)*sin(phi);
M(3,2) = (1-cos(phi))*u(2)*u(3) + u(1)*sin(phi);
M(3,3) = cos(phi) + (1-cos(phi))*u(3)^2;