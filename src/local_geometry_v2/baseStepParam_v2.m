function [shift,slide,rise,tilt,roll,twist,U_mid,r_mid] = baseStepParam_v2(U1,r1,U2,r2,flag)

THRES = 1e-10;

hinge = cross(U1(:,3),U2(:,3));
if(norm(hinge)>THRES)
    [phi,~] = vecAngle(U1(:,3),U2(:,3));
    U1_new = rotateHinge(hinge,phi/2) * U1;
    U2_new = rotateHinge(hinge,-phi/2) * U2;
else
    phi = 0;
    U1_new = U1;
    U2_new = U2;
end

% Calculate the middle frame
U_mid = zeros(3);
U_mid(:,3) = U1_new(:,3);

% For two base-pairs, function bisector is obsolete; it does not work if the twist angle > 180d
switch lower(flag)
    case 'basepair'   % For basePairParam_v2.m
        U_mid(:,1) = bisector(U1_new(:,1),U2_new(:,1));
        U_mid(:,2) = bisector(U1_new(:,2),U2_new(:,2));
    case 'basepairstep'
        [psi,~] = vecAngle(U1_new(:,1),U2_new(:,1));
        if(cross(U1_new(:,1),U2_new(:,1))'*U_mid(:,3) < 0)
            psi = 2*pi - psi;
        end
        U_mid(:,1) = rotateHinge(U_mid(:,3),psi/2) * U1_new(:,1);
        U_mid(:,2) = rotateHinge(U_mid(:,3),psi/2) * U1_new(:,2);
    otherwise
        error('Unknown flag');
end
r_mid = (r1+r2)/2;

% Project the pure translation to the middle frame
tmp = U_mid' * (r2-r1);
shift = tmp(1);
slide = tmp(2);
rise = tmp(3);

% The twist
[~,twist] = vecAngle(U1_new(:,2),U2_new(:,2));
twist = twist * sign(cross(U1_new(:,2),U2_new(:,2))'*U_mid(:,3));

% The roll and tilt
if(norm(hinge)>THRES)
    [phaseAngle,~] = vecAngle(hinge,U_mid(:,2));
    sgn = sign(cross(hinge,U_mid(:,2))'*U_mid(:,3));
    %roll = phi*(180/pi) * abs(cos(phaseAngle)) * sgn;
    %tilt = phi*(180/pi) * abs(sin(phaseAngle)) * sgn;
    roll = phi*(180/pi) * cos(sgn*phaseAngle);
    tilt = phi*(180/pi) * sin(sgn*phaseAngle);
else
    roll = 0;
    tilt = 0;
end