function [Slope,Intercept] = AdaptedPassingBablok (x,y)

% output = 1 => Hypothesis of linearity is NOT rejected
% output = 0 => Hypothesis of linearity is rejected

%% Adapted Passing Bablok regression 
% Written in 2017 by Ulrich Hölker, Freiburg
% Based on Andrea Padoan, PADOVA, ITALY, July 2009
% Link: https://de.mathworks.com/matlabcentral/fileexchange/24894-passing-and-bablok-regression?focused=5144148&tab=function&requestedDomain=true
% adapted to datasets which consist of grouped data (eg grouped
% measurements)
% Main difference: Datapoints coming from one and the same group are
% NOT used for regression parameter estimation

% Input variables x and y have to be matrices (one column = one group)
% In case of unequal sized groups columns will be filled up with NaNs
% During the regression, NaNs will be ignored

m = size(x,2);      % number of colums equals no of groups
p = size(x,1);      % number of rows equals no of elements per group

%% Find NaN     
    n = numel(x);    
    n_nan = sum(sum(isnan(x)));
  %  display ('Number of NaNs: ');
  %  display (n_nan);
    
%% Estimation of beta and alpha, according to Theil regression
    
    % Compute Slopes S_ij^kl between all possible pairs of points
    % BUT ignore pairs of points coming from one and the same group
    count = 1;
    for j=1:m       % Run through all columns, ie groups 
        for k=1:p   % Run through all elements of each group 
                for l=(j*p+1):n
                        Sij(count) = (y(k,j) - y(l)) /(x(k,j) - x(l));
                        count = count + 1;           
                end % of l loop
        end % of k loop
    end % of j loop
    
    % Exclude NaNs and sort Sij
    S = Sij(~isnan(Sij));
    S = sort(S, 'ascend');
    N = length(S);
    
    % Calculate K , that is the number of values of Sij with Sij < -1
    K= 0;
    for i=1:N
        if (S(i) < -1)
            K = K +1;
        end
    end

    % b is estimated by the shifted median of the Sij using K as offset
    % If N is odd
    if rem(N, 2) ~= 0
        b = S ( (N+1)/2+K );
    else % if N is even
        b = 0.5 * ( S (N/2+K) + S (N/2+1+K) );
    end

    Slope = b;
      

%% Estimate slope, intercept and 95 % confidence intervals

    % Variance formula for setting with grouped dataset
    % for derivation see master's thesis
    p_i = 0;
    p_i = sum(~isnan(x));
    sum_pi = sum(p_i);
    temp = 0;
    
    % Compute empirical probabilty of overlapping on x-axis of adjacent
    % groups "left to right"
    qhat = zeros(m,1);
    for k = 1:(m-1)     % for every group
        temp_sum = 0;
        for j = 1:p     % for every element of that group
            temp_count = 0;
            temp_count = sum(x(j,k) <= x(:,k+1));     % Number of elements of group k+1 located right of element (j,k)
            temp_sum = temp_sum + temp_count*(p-temp_count);
        end
        qhat(k) = 2/(p*p*(p-1))*temp_sum;
    end
        
    for i=1:m
        temp = 18*(sum_pi - 1)+4*(sum_pi -1)*(sum_pi -2) -3*(sum_pi - p_i(i))*(p_i(i)-1)*qhat(i) - (p_i(i)-1)*(p_i(i)-2);
    end
    APBVar = 1/36*sum_pi*temp;

    % Define a quantile for the construction of two side confidence interval for
    % beta. 1.96 means 95 % CI

    w = 1.96;
  % C = w * sqrt ( (m*p*(m*p-1)*(2*m*p+5) ) /18 ) ; formula for simple
  % setting (see paper Passing & Bablok)
    C = w * sqrt(APBVar);   % new formula (for derivation see master's thesis)
    M1 = (N - C)/2;
    M1 = round(M1);
    M2 = N-M1+1;

    % Slope confidence intervals
    slope_UB = S(M2+K);
    slope_LB = S(M1+K);
    
    % Convert matrices into vectors and drop NaNs
    x = x(1:m*p);
    y = y(1:m*p);
    x = x(~isnan(x));
    y = y(~isnan(y));

    % Estimation of alpha (=a) and alpha CI
    a = median (y-x*b);
    Intercept = a;
    a_LB = median (y-slope_UB*x);
    a_UB = median (y-slope_LB*x);

        
%% Compute Residual Sum of Squares  

    residuals = y - (x*Slope + Intercept);
    RSS = sum((residuals).^2);

    
%% Statistical test of the assumption of linearity

    % Calculate for every (xi,yi) point the distance from the intercept. These 
    % D values are needed for rank (xi,yi)

    for i =1:length(x)
        D (i) = ( y(i)+(1/b)*x(i)-a )/ sqrt (1+1/(b*b));
    end

    % Create an array with D in the 1st row, x in the 2nd and y in the 3rd
    Dxy = D';
    Dxy (:,2) = x;
    Dxy (:,3) = y;

    % Sort the first column of the array. In this way you rank xy for the
    % distance from the intercept.
    Dxy = sortrows (Dxy, 1);

    % Create a new variable with ranked xy
    ranked_xy = Dxy (:,2:3);

    % Calculate the number of l and L. l denote the number of points (xi, yi)
    % with yi > a+bxi and L the number of point with yi < a+bxi.
    l = 0;
    L = 0;
   
    for i =1:length(x)
       if y (i) > a+b*x(i)
           l = l + 1;
       elseif y (i) < a+b*x(i)
           L = L + 1;
       end
    end

   
    % Give a score r to each point of ranked xy. This rules are specified on
    % Passing Bablok paper (1985). r now is ranked, like ranked_xy.
    for i = 1:length(x)
       if ranked_xy(i,2) > a+b*ranked_xy(i,1)
           r (i) = sqrt(L/l);
       elseif ranked_xy(i,2) < a+b*ranked_xy(i,1)
           r (i) = - sqrt(l/L);
       elseif ranked_xy(i,2) == a+b*ranked_xy(i,1)
           r (i) = 0;
       end
    end   
    
    % Calculate cumulative sum of ranked r.
    cumsum_r = cumsum(r);
    % Calculate max absolute value of cumsum of r.
    max_cumsum_r = max(abs (cumsum_r));

%% Display results of statistical test

% comment the whole section to make in run faster
% if you're not interested in showing those results

    % Calculate limits using critical values of Kolmogorov-Smirnov test.
    p1perc = 1.63*sqrt(L+l);
    p5perc = 1.36*sqrt(L+l);
    p10perc = 1.22*sqrt(L+l);

    display('Test for linear assumption: ');
    display('         H0: x and y are linear relationship');
    display( ' ');

    if max_cumsum_r <= p10perc 
        display ('  p > 0.1, so x and y have a linear relationship. Do not reject H0.');
        output = 1;
        crt_value = p10perc
    end

    if (max_cumsum_r >= p10perc) && (max_cumsum_r < p5perc)
        display ('  p > 0.05, so x and y have a linear relationship. Do not reject H0.');
        output = 1;
        crt_value = p10perc
        crt_value = p5perc
    end

    if (max_cumsum_r >= p5perc) && (max_cumsum_r <= p1perc)
        display ('  0.01 < p < 0.05 and Linear relationship between x and y is rejected.');
        output = 0;
        crt_value = p5perc
        crt_value = p1perc
    end

    if (max_cumsum_r >= p1perc)
        display ('  p < 0.01 and Linear relationship between x and y is rejected.');
        output = 0;
        crt_value = p1perc
    end

%% Plot regression graph

% comment the whole section to make in run faster
% if you're not interested in showing those results

    % Plot point of (xi,yi)
    figure ('Name', 'Passing Bablok Regression fit', 'NumberTitle','off');
    scatter (x,y);
    hold on;

    % Fit the curve from lowest value of x to the largest.
    xhat = sort(x);
    for i = 1:length(x)
        yhat(i) = a+b*xhat(i);
    end
    
    % Create a curve for slope_UB (upper 95% CI) and slope_LB(lower 95% CI)
    for i = 1:length(x)
        yhatLB(i) = a_LB+slope_LB*xhat(i);
        yhatUB(i) = a_UB+slope_UB*xhat(i);
    end
        
    % Define max values for axis
    if max(sort(x)) > max(sort(y))
        axis ([0 max(x) 0 max(x)]);
    else
        axis ([0 max(y) 0 max(y)]);
    end
    
    % Plot the line
    plot (xhat, yhat, 'LineWidth', 2, 'Color', 'b');
    plot (xhat, yhatLB, 'LineWidth', 1, 'Color', 'black', 'LineStyle','--');
    plot (xhat, yhatUB, 'LineWidth', 1, 'Color', 'black', 'LineStyle','--');
    title ('Passing Bablok regression', 'FontSize',18, 'Color', 'b', 'FontName', 'Courier','fontweight','b');
    xlabel ('method B, x values');
    ylabel ('method A, y values');
    
    % Plot bisecting line
    if max(xhat) >= max(yhat) 
        line ([0 max(xhat)], [0 max(xhat)], 'Color', 'r', 'LineStyle', ':' , 'LineWidth', 1);
    else
        line ([0 max(yhat)], [0 max(yhat)], 'Color', 'r', 'LineStyle', ':' , 'LineWidth', 1);
    end

%% Plot scatter plot of ranked residual

% comment the whole section to make in run faster
% if you're not interested in showing those results

    for i = 1:length(x)
        residual(i) = Dxy(i,3) - b*Dxy(i,2)-a;
    end

    figure('Name', 'Ranked residual plot', 'NumberTitle','off');
    hold off
    plot (1:length(x), residual, 'bo', 'MarkerSize',5, 'LineWidth', 2);
    hold on
    line ([ 0 length(x)], [0 0], 'Color', 'r', 'LineWidth', 2);
    title ('Residual Plot', 'FontSize',18, 'Color', 'b', 'FontName', 'Courier','fontweight','b');
    ylabel ('Residuals');
    xlabel ('Rank (xi,yi)');
    makePretty;

%% Plot cumsum_r statistics

% comment the whole section to make in run faster
% if you're not interested in showing those results

    hold off;
    figure('Name', 'cumsum statistic', 'NumberTitle','off');
    % Plot cumsum, added of starting 0 value.
    plot(0:length(x), horzcat([0],cumsum_r), 'LineWidth', 2,'Color','blue');
    hold on;
    % title ('Linearity of tests', 'FontSize',18, 'Color', 'b', 'FontName', 'Courier','fontweight','b');
    ylabel ('cumsum');
    xlabel ('Rank (xi,yi)'); 
    line ([ 0 length(x)], [p5perc p5perc], 'Color', 'g', 'LineWidth', 2);
    line ([ 0 length(x)], [-p5perc -p5perc], 'Color', 'g', 'LineWidth', 2);
    line ([ 0 length(x)], [0 0], 'Color', 'r', 'LineWidth', 2);
    plot(0:length(x), horzcat([0],cumsum_r), 'LineWidth', 2,'Color','blue');
    xlim([0 length(x)]);
    ylim([round(-p5perc)-2 round(p5perc)+2]);
    h = legend('cusum(i)','$\pm h_\gamma + \sqrt(l+L)$');
    set(h,'Interpreter','latex')
    makePretty;


end   % of function


