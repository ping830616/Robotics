function yy = getYY(SetpointMat)
v_BC = [0;0];                   %initial and final velocity
n = size(SetpointMat,1);
%Set the results matrix 4(n-1)*1 size
transMat = zeros(2*(n-1),n);    %transform matrix for Boundary conditions
transMat(1) = 1; transMat(end) = 1;

i=2;
for j = 2:n-1
    transMat(i,j) = 1;
    transMat(i+1,j) =1;
    i=i+2;
end

yy=[transMat*SetpointMat;v_BC*ones(1,3);zeros(2*(n-2),3);];