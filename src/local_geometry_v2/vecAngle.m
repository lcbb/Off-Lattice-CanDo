function [rad,deg] = vecAngle(vec1,vec2)

rad = acos(dot(vec1,vec2)/norm(vec1)/norm(vec2));
deg = rad/pi*180;