% M2-EST-TP2
% TP Estimation
% Cree le 25/02/20

close all; % Clear figs
clc; % Clear term

disp('M2-EST-TP2');

%display_scenarios;
nums = [1, 3, 4, 5];
% Note : pour scenar 4, f1 et f4
% Note : pour scenar 5, f1, f2, f3, f4
for num = nums
    filepath = sprintf("scenario/scenario%dbis.mat", num);
    load(filepath);

    y = [yobs1, yobs2];
    ypt = [yobs1d, yobs2d];
    [F1, F2] = moindres_carres(y, ypt, u, 0);
    fprintf("\n Scenario %d :\n", num);
    
    fprintf('Eq%d f%d = %.2f\n', 1, 1, F1(3));
    fprintf('Eq%d f%d = %.2f\n', 1, 2, F1(2));
    fprintf('Eq%d f%d = %.2f\n', 1, 4, F1(1));
    
    fprintf('Eq%d f%d = %.2f\n', 2, 2, F2(1));
    fprintf('Eq%d f%d = %.2f\n', 2, 3, F2(2));
end

% Q3
% On constate que les courbes très différentes de zéros sont celles qui ont
% des défauts. On peut donc dire que le scénario 1 est celui sans défaut et
% le scénario 4 qui a un défaut sur r1 et les deux autres on des défaut sur
% r1 et r2.

function [F1, F2] = moindres_carres(y, ypt, u, ~)
    nb_vals = size(y, 1);
    % AX = B => X = pinv(A) * B
    
    % --- Eq1 ---
    % Error : here
    A1 = zeros(nb_vals, 3);
    A1(:,1) = 1;                                         % 0.3 * f4 (f4 - 2)(y1 - u - f1 - f2)
    A1(:,2) = - 2 * ypt(:,1) - 0.3;                      % f2
    A1(:,3) = -0.3;                                      % f1
    B1 = - (-0.3 * u + (0.3 + 2 * ypt(:,1)) .* y(:,1));
    F1 = pinv(A1) * B1;
    
    f1 = F1(3);
    f2 = F1(2);
    rs = roots([-0.3 * (f1 + f2), 0.6 * (f1 + f2), -F1(1)]); 
    f4 = NaN;
    for r = rs'
        if (0 <= r) && (r <= 1)
            f4 = r(1);
        end
    end
    %fprintf('Roots = %s\n', mat2str(rs));
    F1(1) = f4;
    
    % --- Eq2 ---
    A2 = zeros(nb_vals, 2);
    A2(:,1) = 0.3;                                       % f2
    A2(:,2) = -2 * ypt(:,2) - 0.3;                       % f3
    B2 = - (2 * ypt(:,2) .* y(:,2) - 0.3 * y(:,1) + 0.3 * y(:,2));
    F2 = pinv(A2) * B2;
end 

function display_scenarios
    nums = [1, 3, 4, 5];
    i = 1;
    for num = nums
        filepath = sprintf("scenario/scenario%dbis.mat", num);

        load(filepath);

        [eq1, eq2] = rra_sd([yobs1, yobs2], [yobs1d, yobs2d], u, 0);
        subplot(4, 2, i);
        plot(1:size(eq1,1), eq1);
        title(sprintf('Mean eq%d : %.5f', 1, mean(eq1)));

        subplot(4, 2, i+1);
        plot(1:size(eq2,1), eq2);
        title(sprintf('Mean eq%d : %.5f', 2, mean(eq2)));

        i = i + 2;
    end
end

function [eq1_sd, eq2_sd] = rra_sd(y, ypt, u, ~)
    eq1_sd = -0.3 * u + (0.3 + 2 * ypt(:,1)) .* y(:,1);
    eq2_sd = 2 * ypt(:,2) .* y(:,2) - 0.3 * y(:,1) + 0.3 * y(:,2);
end