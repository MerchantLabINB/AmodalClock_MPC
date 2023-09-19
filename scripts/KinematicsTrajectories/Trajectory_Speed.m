function Speed = Trajectory_Speed(Trajectory)
% Calculation of Trajectory speed in the first 3 principal component
% Input:
%       Trajectory -> Matrix -> NDim X NBin
%           Trajectory Time Series        
% Outút:
%       Speed -> Matrix -> NBin X 1
% Time differential (defaul value 1)
dt = .20; %1  .020
%Only the 3 principal component are taken
Trajectory = Trajectory(1:3, :)';
% Calculation of velocity (Vector)
Velocity = diff(Trajectory) / dt;
% Calculation of speed
Speed = sqrt(sum(Velocity.^2, 2));


% % pc  =2;
% scale = 0.4;
% x = 1:50;
% y = Speed(x,1)';
% % dy = Velocity(x,pc)';
% % dx = gradient(x);
% figure; hold on;
% plot(x,y)
% quiver(x,y,x,-y,scale,'w','LineWidth',1)
% % hold on; plot( x, y)