% This is for planning the cubic spline trajectory for a
% RRR planar robotic arm
% Inputs:
% - Goal matrix             G (via points with designated time of arrival)
% - Cartesian or Joint?     CJ
% - Link Info               l1 l2 l3

%% Settings
addpath(genpath(pwd));  %add subfolder path to search path

% Link length
l1 = 5;
l2 = 3;
l3 = 1;

% CJ: 0: in Cartesian space
%     1: in joint space
CJ = 0;

%       t    x    y    theta
G = [   0   -4    0      120 ;
        2   -5    5      45 ;
        4    2    3      30 ;
        9    2   -3       0 ];

% End of Settings
%% Animation and video recording settings

% Animation?
Ani = 1;
%video recording
videoSave = 0;
videoSave = videoSave*Ani;
videoName = 'cubic-cart';
if videoSave == 1
    AVI1 =VideoWriter(videoName,'MPEG-4');
    open(AVI1);
end

%%

%Line and Point width
OW = 2;
LW = 4;
%total number of steps
step = size(G,1);
dt = diff(G(:,1));
Tint = 0.05;
N = dt/Tint+1;

% 由matrix抓資料出來，單純只是讓後續code比較易讀
To = G(:,1);
Q(:,1) = G(:,2);
Q(:,2) = G(:,3);
Q(:,3) = G(:,4);


%% To solve the coefficents all cubic splines
% DOF = T * A
% function getYY: to get the vector DOF
if CJ == 0
    yy = getYY(G(:,2:end)); 
    XG = yy(:,1);
    YG = yy(:,2);
    DG = yy(:,3);
elseif CJ == 1
    q = InverseKin(G(:,2),G(:,3),G(:,4),l1,l2);
    q1 = q(:,1);
    q2 = q(:,2);
    q3 = q(:,3);
    yy = getYY([q1 q2 q3]);
    XG = yy(:,1);
    YG = yy(:,2);
    DG = yy(:,3);
else
    disp('Wrong index!');
end

% function genCubicMatrix: to generate the massive Matrix T, 
% for solving all coeff. of cubic equations
M = genCubicMatrix(G);

XA = M\XG; %inv(M)*XG;
YA = M\YG; %inv(M)*YG;
DA = M\DG; %inv(M)*DG;

X = [];
Y = [];
D = [];
T = [];
for i = 1:step-1
    t = linspace(0,dt(i),N(i))';
    if i==1
        T = t;
    else
        T = [T;t+T(end)];
    end
    %obtain every point of X Y D
    X = [X; XA(4*(i-1)+1,1)*ones(N(i),1)+XA(4*(i-1)+2,1)*t+XA(4*(i-1)+3,1)*t.^2+XA(4*(i-1)+4,1)*t.^3];
    Y = [Y; YA(4*(i-1)+1,1)*ones(N(i),1)+YA(4*(i-1)+2,1)*t+YA(4*(i-1)+3,1)*t.^2+YA(4*(i-1)+4,1)*t.^3];
    D = [D; DA(4*(i-1)+1,1)*ones(N(i),1)+DA(4*(i-1)+2,1)*t+DA(4*(i-1)+3,1)*t.^2+DA(4*(i-1)+4,1)*t.^3];
end

%% Plot graphs
if CJ == 0
    figure(1);
    subplot(1,3,1);
    plot(T,X,'LineWidth',LW);xlabel('t');ylabel('X(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),G(:,2),'ro','LineWidth',OW);hold off;
    subplot(1,3,2);
    plot(T,Y,'LineWidth',LW);xlabel('t');ylabel('Y(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),G(:,3),'ro','LineWidth',OW); hold off;
    subplot(1,3,3);
    plot(T,D,'LineWidth',LW);xlabel('t');ylabel('D(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),G(:,4),'ro','LineWidth',OW); hold off;
    
elseif CJ == 1
    
    %obtain the links' X, Y, D
    CX = l1*cos(X)+l2*cos(X+Y);
    CY = l1*sin(X)+l2*sin(X+Y);
    CD = (X + Y+D).*(180/pi);
    
    figure(11);
    subplot(1,3,1);
    plot( T,X.*(180/pi),'LineWidth',LW );xlabel('t');ylabel('q1(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),q1.*(180/pi),'ro','LineWidth',OW);
    hold off;
    subplot(1,3,2);
    plot( T,Y.*(180/pi),'LineWidth',LW );xlabel('t');ylabel('q2(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),q2.*(180/pi),'ro','LineWidth',OW); hold off;
    subplot(1,3,3);
    plot(T,D.*(180/pi),'LineWidth',LW );xlabel('t');ylabel('q3(t)');set(gca,'fontsize',25);
    hold on;plot(G(:,1),q3.*(180/pi),'ro','LineWidth',OW);
    hold off;
    
end

if CJ == 0
    figure(2);
    plot(X,Y,'LineWidth',LW);axis equal;grid on;grid minor;set(gca,'fontsize',25);
    xlim([-(l1+3) (l1+3)]);ylim([-(l1+3) (l1+3)]);
    xlabel('X');ylabel('Y');
    hold on;plot(G(:,2),G(:,3),'ro','LineWidth',OW); hold off;
elseif CJ == 1
    figure(12);
    subplot(1,2,1)
    plot(CX,CY,'LineWidth',LW);axis equal;grid on;set(gca,'fontsize',25);
    xlim([-(l1+3) (l1+3)]);ylim([-(l1+3) (l1+3)]);
    xlabel('X');ylabel('Y');
    hold on;plot(G(:,2),G(:,3),'ro','LineWidth',OW); hold off;
    subplot(1,2,2);
    plot(T, CD,'LineWidth',LW);
    xlabel('t');ylabel('phi(t)');
    grid on;hold on;plot(G(:,1),G(:,4),'ro','LineWidth',OW); hold off;
    
end

%% Animation
N = length(X);
if Ani == 1
    endP = N;
else
    endP = 1;
end
if CJ == 0
    fig = figure(3);
    grid on;grid minor;
    q = InverseKin(X,Y,D,l1,l2);
    q1 = q(:,1);
    q2 = q(:,2);
    q3 = q(:,3);
    QQ(:,1) = X;
    QQ(:,2) = Y;
elseif CJ == 1
    QQ(:,1) = CX;
    QQ(:,2) = CY;
    q1 = X;
    q2 = Y;
    q3 = D;
    fig = figure(13);
    set(13, 'Position', [300 200 560 420]);
end

grid on;grid minor;
plot(QQ(:,1),QQ(:,2),'LineWidth',LW);             %plot the trajectory
set(gca,'fontsize',25);
axis equal;grid on;grid minor;
xlim([-(l1+l2) (l1+l2)]);ylim([-(l1+l2) (l1+l2)]);
xlabel('X');ylabel('Y');
hold on;
plot(Q(:,1),Q(:,2),'r.','MarkerSize',25);	% 畫出initial, via, and final points
d = linspace(0,2*pi);
plot( (l1-l2)*cos(d), (l1-l2)*sin(d), 'm');	% 畫出work space，inner boundary
plot( (l1+l2)*cos(d), (l1+l2)*sin(d), 'm');	% 畫出work space，outer boundary
plot(0,0,'ko');                             % 畫出圓當成地桿
tic;
XX = [ zeros(N,1) l1*cos(q1) l1*cos(q1)+l2*cos(q1+q2) l1*cos(q1)+l2*cos(q1+q2)+l3*cos(q1+q2+q3) ];
YY = [ zeros(N,1) l1*sin(q1) l1*sin(q1)+l2*sin(q1+q2) l1*sin(q1)+l2*sin(q1+q2)+l3*sin(q1+q2+q3) ];
for i = 1:endP
    f1 = plot(XX(i,:),YY(i,:),'k','linewidth',4);   % 畫出手臂
    pause(0.08);
    if (i<N) && (Ani == 1)
        delete(f1);
    end
    if videoSave == 1
        frame = getframe(fig);
        writeVideo(AVI1,frame);
    end
end
hold off;

if videoSave==1
    close(AVI1);
end