function [coef,resid,ypred] = curvefitting6functions2019(x,y,model,task)
%%HM 3/18/2015

DFE = numel(x)-2;
DFR = 1;
warning off;

y(y(:) == 0) = 0.000001;

switch(model)
    case(1) %linear
        %Linear Curve Fit: y = mx+b
        p = polyfit(x, y, 1);
        m = p(1);
        b = p(2);
        J = sum((m*x+b-y).^2); %sse
        mu = mean(y);
        S = sum((y-mu).^2); %sst
        SSR = S-J; %SS regression
        r2 = 1-J/S;
        mse = J/DFE;
        msr = SSR/DFR;
        F = msr/mse;                                 %	Estadistico	F
        P = 1 - fcdf(F, DFR, DFE);                   %   Estadistico P
        
        ypred = (m*x+b);  %predicted
        resid = y - ypred;
        
    case(2) %power
        %Power Function Curve Fit: y = bx^m
        p = polyfit(log10(x), log10(y), 1);
        m = p(1); %slope
        b = 10^p(2); %cte
        J = sum((b.*x.^m-y).^2);
        mu = mean(y);
        coef(1) = m;
        coef(2) = b;
        S = sum((y-mu).^2);   %SST
        
        ypred = (b.*x.^m);%polyval(coef,log10(x));  %predicted
        resid = y - ypred;
        
        SSR = sum(abs(log10(ypred)-mean(log10(y))).^2);   %Sum of squares regression
        SSE = sum(abs(log10(ypred)-log10(y)).^2);         %Sum of squares error
        
        r2 = SSR/(SSE+SSR);
        mse = SSE/DFE;
        msr = SSR/DFR;
        F = msr/mse;                                 %	Estadistico	F
        P = 1 - fcdf(F, DFR, DFE);                   %   Estadistico P
    case(3) %expontial
        %Exponential Curve Fit: y = b^mx
        p = polyfit(x, log10(y), 1);
        m = p(1); %slopw
        b = 10^p(2);
        J = sum((b.*10.^(m*x)-y).^2);
        mu = mean(y);
        S = sum((y-mu).^2);
        
        ypred = (b.*10.^(m*x));
        resid = y - ypred;
        
        SSR = sum(abs(log10(ypred)-mean(log10(y))).^2);   %Sum of squares regression
        SSE = sum(abs(log10(ypred)-log10(y)).^2);         %Sum of squares error
        
        r2 = SSR/(SSE+SSR);
        
        mse = SSE/DFE;
        msr = SSR/DFR;
        F = msr/mse;
        P = 1 - fcdf(F, DFR, DFE);                   %   Estadistico P
    case(4) %logarithmic
        %Logatithmic Currve Fit: b + mlnx
        xlog = log(x);
        p = polyfit(xlog, y, 1);
        m = p(1);
        b = p(2);
        ypred = m*xlog + b;
        resid = y - ypred;
        J = sum((y-ypred).^2);
        mu = mean(y);
        S = sum((y-mu).^2);
        
        
        SSR = sum((ypred-mean(y)).^2);   %Sum of squares regression
        SSE = sum((ypred-y).^2);         %Sum of squares error
        
        r2 = SSR/(SSE+SSR);
        
        mse = SSE/DFE;
        msr = SSR/DFR;
        F = msr/mse;
        P = 1 - fcdf(F, DFR, DFE);                   %   Estadistico P
        
    case (5)%%gaussian
        
        gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
        plotfigure = 1;
       
       % [mx,my,sdy] = plotRegression(x,y,task,plotfigure);        
        [maxyy,inxy] = max(y);
       % initmean = mx(inxy);
        initmean = mean(x);
        if task < 3
           initsd = x(end)/4;%(mx(4) - mx(1));
        else
            initsd = x(end)/6;%(mx(4) - mx(1));
        end
        initheight = maxyy+2;%/2;
        initcte = min(y)-2;
        if initcte<0
            initcte = 0;
        end
        BoundWidht = 0.8;
        BoundWidht2 = 0.8;
        BoundWidht3 = 0.2;
        Boundheight = 1.8;
        
        lowermean = x(1) - x(1)*BoundWidht;
        uppermean = initmean+initmean*BoundWidht;
        lowersd = initsd - initsd*BoundWidht2;
        uppersd = initsd + initsd*BoundWidht;
        lowerheight = -1;
        upperheight = initheight + initheight*Boundheight;
        lowercte = 0;
        uppercte = initcte + initcte*BoundWidht3;
        
        
        
        
       % fopts = fitoptions('Method','LinearLeastSquares', 'Normalize','on','Robust','on', 'Lower', [lowerheight lowermean lowersd lowercte], 'Upper', [upperheight uppermean uppersd uppercte]);
        %StartPoints = [initheight initmean initsd initcte];
        
     %   [fg1,gofg1,out1] = fit(x,y,'gauss1',fopts);
      %  [fg,gofg] = fit(x',y',gaussEqn,'Start',[initheight initmean initsd initcte],'Lower', [lowerheight lowermean lowersd lowercte],'Upper', [upperheight uppermean uppersd uppercte]);
        [fg,gofg] = fit(x,y,gaussEqn,'Start',[initheight initmean initsd initcte],'Lower', [lowerheight lowermean lowersd lowercte],'Upper', [upperheight uppermean uppersd uppercte]);
        
        ypred = feval(fg,x);
     %   plot(x,ypred);
%        hold off;
        resid = ypred-y;
        b = fg.b; %mean
        m = fg.c; %sd
        r2 = gofg.adjrsquare;
        P = gofg.rsquare;
        c = true;
        coef(1,5) =  fg.a+ fg.d; %height
   
        
    case (6)
        
%         plot_flag = 0;
%         miny = min(y);
%         maxy = max(y);
%         x50 = mean(x);
%         slope = 1;
%         fixed_params=[miny, maxy , x50 , [NaN]];
%         initial_params =[miny, maxy , x50 , slope];
%         [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag);
%          ypred = stat.ypred;
%          c = true;
%          b = stat.param;
%          m = stat.paramCI(1);
%          P = stat.paramCI(2);
         
         
        [ypred,p] = fit_logistic(x,y);
        b = p(1);
        m = p(2);
        P = p(3);
        c = true;
        resid = y - ypred;        
        [r2 rmse] = rsquare(y,ypred,c);
        
        %         numrep = max(y);
        %         totrep = ones(numel(x),1)*numrep;
        %         [bT,dev,stats] = glmfit(x,[y totrep],'binomial');
        %         ypred = glmval(bT,x,'logit');
        %         ypred = ypred*numrep;
        %         resid = ypred-y;
        %         b = bT(1);
        %         m = bT(2);
        %         c = true;
        %         %r2 = stats.s;
        %         P = stats.p(2);
        %        [r2 rmse] = rsquare(y,ypred,c);
        
        %         %.21
        %         x21(1,1) = (log(0.26588)-b(1,1))/b(2,1);
        %         %.79
        %         x79(1,1) = (log(3.76278)-b(1,1))/b(2,1);
        %         %.5
        %         x50(1,1) = (log(1)-b(1,1))/b(2,1);
        
    otherwise
        disp('Unknown method.')
end

coef(1,1) = b;
coef(1,2) = m;
coef(1,3) = r2;
coef(1,4) = P;

return 

% X = [ones(size(x)) log10(x)]; %regression matrix
% [b,bint,r,rint,stats] = regress(log(y),X,0.95);




function [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
% Optimization of parameters of the sigmoid function
%
% Syntax:
%       [param]=sigm_fit(x,y)       
%
%       that is the same that
%       [param]=sigm_fit(x,y,[],[],[])     % no fixed_params, automatic initial_params
%
%       [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
%       [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor
%       [param]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
% param = [min, max, x50, slope]
%
% if fixed_params=[NaN, NaN , NaN , NaN]        % or fixed_params=[]
% optimization of "min", "max", "x50" and "slope" (default)
%
% if fixed_params=[0, 1 , NaN , NaN]
% optimization of x50 and slope of a sigmoid of ranging from 0 to 1
%
%
% Additional information in the second output, STAT
% [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
%
% Example:
% %% generate data vectors (x and y)
% fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
% param=[0 1 5 1];  % "min", "max", "x50", "slope"
% x=0:0.1:10;
% y=fsigm(param,x) + 0.1*randn(size(x));
%
% %% standard parameter estimation
% [estimated_params]=sigm_fit(x,y)
%
% %% parameter estimation with forced 0.5 fixed min
% [estimated_params]=sigm_fit(x,y,[0.5 NaN NaN NaN])
%
% %% parameter estimation without plotting
% [estimated_params]=sigm_fit(x,y,[],[],0)
%
%
% Doubts, bugs: rpavao@gmail.com
% Downloaded from http://www.mathworks.com/matlabcentral/fileexchange/42641-sigmoid-logistic-curve-fit

% warning off

x=x(:);
y=y(:);

if nargin<=1 %fail
    fprintf('');
    help sigm_fit
    return
end

automatic_initial_params=[quantile(y,0.05) quantile(y,0.95) NaN 1];
if sum(y==quantile(y,0.5))==0
    temp=x(y==quantile(y(2:end),0.5));    
else
    temp=x(y==quantile(y,0.5));
end
automatic_initial_params(3)=temp(1);

if nargin==2 %simplest valid input
    fixed_params=NaN(1,4);
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==3
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==4
    plot_flag=1;    
end

if exist('fixed_params','var')
    if isempty(fixed_params)
        fixed_params=NaN(1,4);
    end
end
if exist('initial_params','var')
    if isempty(initial_params)
        initial_params=automatic_initial_params;
    end
end
if exist('plot_flag','var')
    if isempty(plot_flag)
        plot_flag=1;
    end
end

%p(1)=min; p(2)=max-min; p(3)=x50; p(4)=slope como em Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
%f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));

f_str='f = @(param,xval)';
free_param_count=0;
bool_vec=NaN(1,4);
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        f_str=[f_str ' param(' num2str(free_param_count) ')'];
        bool_vec(i)=1;
    else
        f_str=[f_str ' ' num2str(fixed_params(i))];
        bool_vec(i)=0;
    end
    if i==1; f_str=[f_str ' + (']; end
    if i==2;
        if isnan(fixed_params(1))            
            f_str=[f_str '-param(1) )./ (   1 + 10.^( (']; 
        else
            f_str=[f_str '-' num2str(fixed_params(1)) ')./ (1 + 10.^((']; 
        end
    end    
    if i==3; f_str=[f_str ' - xval ) *']; end
    if i==4; f_str=[f_str ' )   );']; end
end

eval(f_str)

[BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,initial_params(bool_vec==1));
stat.param=BETA';

% confidence interval of the parameters
stat.paramCI = nlparci(BETA,RESID,'Jacobian',J);

% confidence interval of the estimation
[stat.ypred,delta] = nlpredci(f,x,BETA,RESID,'Covar',COVB);
stat.ypredlowerCI = stat.ypred - delta;
stat.ypredupperCI = stat.ypred + delta;

% plot(x,y,'ko') % observed data
% hold on
% plot(x,ypred,'k','LineWidth',2)
% plot(x,[lower,upper],'r--','LineWidth',1.5)

free_param_count=0;
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        param(i)=BETA(free_param_count);
    else
        param(i)=fixed_params(i);
    end    
end
    
if plot_flag==1 
    x_vector=min(x):(max(x)-min(x))/100:max(x);
    plot(x,y,'k.',x_vector,f(param(isnan(fixed_params)),x_vector),'r-')
    xlim([min(x) max(x)])
end

return

