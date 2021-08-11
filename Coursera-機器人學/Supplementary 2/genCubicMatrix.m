%this is to construct the T matrix to solve the cubic equation coefff for
%via point n >= 1
%T matrix = 4(n-1)*4(n-1) size
function  T_Mat = genCubicMatrix(G)

% Example of G
% G = [    0   -4    0      90 ;
%          2    0    3      45 ;
%          4    3    3      30 ;
%          7    4    0       0 ];
%Setpt_X = G(:,2)';           % n elements
%Setpt_Y = G(:,3)';            % n elements
%Setpt_theta = G(:,4)';     % n elements
TimeEle = G(:,1)';            % n elements also

% Setpt_X = [-4 0 3 4];           % n elements
% Setpt_Y = [0 3 3 0];            % n elements
% Setpt_theta = [90 45 30 0];     % n elements
% TimeEle = [0 2 4 7];            % n elements also


n = numel(TimeEle);
dt_size = numel(TimeEle);
dt = zeros(n-1,1);
for i = 1:dt_size-1
    dt(i) = TimeEle(i+1)-TimeEle(i);    %matrix of each time inteval obtained
end



%result matrix with boundaries, special condition*2, 
%velocity and acc continuity 
%---------------------------------------------------------------
%Set the T matrix 4(n-1)*4(n-1) size
%3(n-1)
T_Mat = zeros(4*(n-1),4*(n-1));
i=1;j=1;
for k=1:n-1
    T_Mat(i,j) = 1;
    T_Mat(i+1,j:j+3) = [1 dt(k) dt(k)^2 dt(k)^3];
    i = i+2;
    j = j+4;
end
i=1+2*(n-1);j=1;

T_Mat(i,2) = 1;
T_Mat(i+1,end-2:end) = [1 2*dt(end) 3*dt(end)^2];
i = i+2;
for k=1:n-2
    T_Mat(i,j:j+7) = [0 1 2*dt(k) 3*dt(k)^2 0 -1 0 0];
    T_Mat(i+1,j:j+7) = [0 0 2 6*dt(k) 0 0 -2 0;];
    i = i+2;
    j = j+4;
end

end
