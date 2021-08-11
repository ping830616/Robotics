%Inverse Kinematics particularly for RRR planer manipulator
%Inputs: X, Y, Theta(in degrees), l1 and l2 length
function q = InverseKin(Xcoor, Ycoor, theta, l1, l2)

n = size(Xcoor,1);
q = zeros(n,3);
q(:,2) = acos((Xcoor.^2+Ycoor.^2-l1^2*ones(n,1)-l2^2*ones(n,1))./(2*l1*l2));
psi = acos((Xcoor.^2+Ycoor.^2+l1^2*ones(n,1)-l2^2*ones(n,1))./(2*l1*sqrt(Xcoor.^2+Ycoor.^2)));
if q(:,2) < 0
    q(:,1) = atan2(Ycoor,Xcoor)+psi;
elseif q(:,2) >0
    q(:,1) = atan2(Ycoor,Xcoor)-psi;
end
q(:,3) = theta.*(pi/180)-q(:,1)-q(:,2);
