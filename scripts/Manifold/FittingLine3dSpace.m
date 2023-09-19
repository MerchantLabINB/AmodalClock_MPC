
%Duration Projection
%x = TapsPry(:,1);
%y = TapsPry(:,2);
%z = TapsPry(:,3);



function [xfit,yfit,zfit,Vp] = FittingLine3dSpace(x,y,z)

%Taps Duration Plane
% x %Manifold.TapsReps(:,1);
% y %Manifold.TapsReps(:,2);
% z %Manifold.TapsReps(:,3);

% xyz - all Taps
% [V, D2] = eig(cov(xyz));
xyz=[x y z];
r=mean(xyz);
xyz=bsxfun(@minus,xyz,r);

%% Method 1
% [V, D] = eig(cov(xyz));
% Vp = V(:,1);
% Var_Axes_Trj = diag(D) / sum(diag(D));

%% Method 2
[~,D,V]=svd(xyz,0);
Vp = V(:,3);
Var_Axes_Trj = diag(D) / sum(diag(D));


% Var_Axes_Trj = diag(D) / sum(diag(D));
%Line projection
%Method 1
% x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(1,1)*V(1,1);
% y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(1,1)*V(2,1);


%Method 2
x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(3,1)*V(1,1);
y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(3,1)*V(2,1);

xfit = x_fit(z);
yfit = y_fit(z);
zfit = z;
hold on
plot3(xfit,yfit,z,'Color',[.5,.5,.5],'linewidth',3)%6 -Pry2-3
% ylabel('Pry-2 (a.u.)')
% zlabel('Pry-3 (a.u.)')
% hold on;
% plot3(xfit(1),yfit(1),z(1),'ro','MarkerSize',30)%6 -Pry2-3

%PCA
% view([-90,0])
% ylim([-60 180])%Pry-2
% zlim([-30 200])%Pry-3
% set(gca, 'ZTick', 0:50:200)
%-----GPFA-----%
% view([90,0])

% figure(1),clf(1)
% hold on;
% plot3(x,y,z,'.w')
% plot3(x_fit(z),y_fit(z),z,'r','linewidth',2)
%Duration Plane
% [V, D]=eig(cov(xyz));
% NOTES:
% C=(xyz'*xyz)/(500);           % variance-covariance matrix of X
% [R,D]=svd(C,0);
% 1) Direction of best fit line corresponds to R(:,1)
% 2) R(:,1) is the direction of maximum variances of dX
% 3) D(1,1) is the variance of dX after projection on R(:,1)
% 4) Parametric equation of best fit line: L(t)=X_ave+t*R(:,1)', where t is a real number
% 5) Total variance of X = trace(D)
% Coefficient of determineation; R^2 = (explained variance)/(total variance)
% D=diag(D);
% R2=D(1)/sum(D);

% https://www.codefull.net/2015/06/3d-line-fitting/