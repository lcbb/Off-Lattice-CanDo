function umid = bisector(u1,u2)

hinge = cross(u1,u2);
[phi,~] = vecAngle(u1,u2);

umid = rotateHinge(hinge,phi/2) * u1;