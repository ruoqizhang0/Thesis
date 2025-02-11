clear, clc, close all

z = chebpts(512);
    FILES = dir("/Users/alagris/Documents/MATLAB/results/*.mat");
    SIZE = size(FILES);
    LB = zeros(SIZE);
    Pe = zeros(SIZE);
    A = zeros(SIZE);
    B = zeros(SIZE);
    K = zeros(SIZE);
    T = zeros(SIZE);
    k_critical = cell(length(FILES), 1);
    %k_critical = zeros(SIZE);
    for i = 1:size(FILES,1)
        x = FILES(i).name;
        %taget = {'results_Pe_1.000000e-01.mat',...
        %    'results_Pe_8.697490e-01.mat'};
        %if ismember(x, taget)
           disp(x);
           load(fullfile('/Users/alagris/Documents/MATLAB/results', x)) 
           k_critical{i} = find(output_hr.pres<=1e-5);
           LB(i) = output_hr.LB;
           Pe(i) = output_hr.Pe;
           A(i) = output_hr.a;
           B(i) = output_hr.b;
           X{i} = output_hr.xi;
           T(i) = output_hr.sol.solvertime;
           xi(:,i) = chebpolyval(flipud(X{i}), z);
        %end
    end
    %Pe = Pe(Pe ~= 0);
    %X = X(Pe ~= 0);
    %xi = xi(:, Pe ~= 0);
    %sort
    [Pe, sorting] = sort(Pe);
    A = A(sorting);
    B = B(sorting);
    LB = LB(sorting);
    T = T(sorting);
    k_critical = k_critical(sorting);
    X = X{sorting};

    figure();
    %plot(z, xi, '.', 'MarkerSize', 10);
    for i = 1:length(Pe)
        loglog(Pe(i), k_critical{i}, '.', 'MarkerSize', 10);
        hold on
    end
    %loglog(Pe, k_critical, '.', 'MarkerSize', 10);
    %hold on
    %loglog(Pe, LB.*Pe, '-.', 'LineWidth', 4, 'Color', 'r');
    %loglog(Pe, LB.*Pe, '-.', 'MarkerSize', 20);
    %xlabel('$\it{z}$', 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('$\it{Pe}$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 18);
    %title('Lower bound*$\it{Pe}$ vs $\it{Pe}$', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 1.5);
    box on;