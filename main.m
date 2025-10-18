%This code repository can only be used for related research and cannot be used for commercial purposes.
%For related matters, please contact cwd2018@hhu.edu.cn
clc;
clear;
close all;

%% parameter
SearchAgents_no = 30;    % population size
Max_iteration = 500;     % Maximum iterations
num_runs = 30;           % Number of runs per algorithm per function
Function_name=1;         % FunctionName of cec2020 link Get_Functions_cec2020
%Function_name="F1";     % FunctionName of 23 test functions link Get_Functions_details
%Function_name="three_bar";   % FunctionName of Engineering problems link Get_Functions_details
dim=10;
[lb,ub,dim,fobj]=Get_Functions_cec2020(Function_name,dim); % Getting the test function if cec2020
%[lb,ub,dim,fobj]=Get_Functions_details(Function_name);  % Getting the test function else
lb = lb .* ones(1, dim);
ub = ub .* ones(1, dim);

% Define the list of algorithms
algorithms = {
    @MSFOA, "MSFOA";
    @SFOA, "SFOA";
    @HGWOSCA, "HGWOSCA";
    @AOA, "AOA";
    @PLO, "PLO";
    @RIME, "RIME";
    @SFOA, "CFOA";
};

num_algorithms = size(algorithms, 1);      

% Initialization result storage
results = struct();
for i = 1:num_algorithms
    results(i).name = algorithms{i,2};
    if strcmp(results(i).name, 'PLO')
        results(i).all_curves = cell(num_runs, 1); % Store all running curves (cell array)
    else
        results(i).all_curves = zeros(num_runs, Max_iteration); % Store all running curves (matrix)
    end
end

for algo_idx = 1:num_algorithms
    algo_func = algorithms{algo_idx,1};
    algo_name = algorithms{algo_idx,2};
    
    for run = 1:num_runs
        try
            if strcmp(algo_name, 'MSFOA')
                [~, ~, Best_curve] = algo_func(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);
            elseif strcmp(algo_name, 'SFOA')
                [~, ~, Best_curve] = algo_func(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);
            elseif strcmp(algo_name, 'HGWOSCA')
                [~, ~, Best_curve] = algo_func(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
            elseif strcmp(algo_name, 'AOA')
                [~, ~, Best_curve] = algo_func(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
            elseif strcmp(algo_name, 'PLO')
                [~, ~, Best_curve] = algo_func(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
                results(algo_idx).all_curves{run} = Best_curve;
                continue; % Skip subsequent store operations
            elseif strcmp(algo_name, 'RIME')
                [~, ~, Best_curve] = algo_func(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
            elseif strcmp(algo_name, 'CFOA')
                [~, ~, Best_curve] = algo_func(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
            end
            results(algo_idx).all_curves(run,:) = Best_curve(1:min(end, Max_iteration));
        catch ME
            disp([ME.message]);
            if strcmp(algo_name, 'PLO')
                results(algo_idx).all_curves{run} = []; % Storing an Empty Array
            else
                results(algo_idx).all_curves(run,:) = Inf;  % If something goes wrong, set to infinity
            end
        end
    end
end

% Calculate the global maximum number of iterations
max_global_iter = Max_iteration;
for run = 1:num_runs
    plo_curve = results(5).all_curves{run};
    if ~isempty(plo_curve)
        max_global_iter = max(max_global_iter, length(plo_curve));
    end
end

% Creating a Graphics Window
figure('Units', 'centimeters', 'Position', [0 0 16 12], 'Color', 'w');

hold on;
colors = lines(num_algorithms); 
line_styles = {'-', '--', ':', '-.', '-', '--', ':'};

% Processing and plotting curves for each algorithm
for algo_idx = 1:num_algorithms
    algo_name = results(algo_idx).name;
    padded_curves = zeros(num_runs, max_global_iter);
    
    if strcmp(algo_name, 'PLO')
        for run = 1:num_runs
            curve = results(algo_idx).all_curves{run};
            if isempty(curve)
                padded_curves(run, :) = NaN; % Flagging invalid data
            else
                curve_len = length(curve);
                % Dynamic Fill Curve
                if curve_len < max_global_iter
                    padded_curves(run, 1:curve_len) = curve;
                    padded_curves(run, curve_len+1:end) = curve(end);
                else
                    padded_curves(run, :) = curve(1:max_global_iter);
                end
            end
        end
    else
        for run = 1:num_runs
            curve = results(algo_idx).all_curves(run, :);
            curve_len = length(curve);
            
            if curve_len < max_global_iter
                padded_curves(run, 1:curve_len) = curve;
                padded_curves(run, curve_len+1:end) = curve(end);
            else
                padded_curves(run, :) = curve(1:max_global_iter);
            end
        end
    end
    
    % Calculate the average curve (ignoring NaN values)
    avg_curve = mean(padded_curves, 1, 'omitnan');
    plot(1:max_global_iter, avg_curve, ...
        'LineWidth', 1.5, ...
        'Color', colors(algo_idx,:), ...
        'LineStyle', line_styles{algo_idx}, ...
        'DisplayName', results(algo_idx).name);
end

xlabel('Iteration', 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Fitness Value', 'FontSize', 8, 'FontWeight', 'bold');
title(['Algorithm Comparison - F' num2str(Function_name)], ...
      'FontSize', 9, 'FontWeight', 'bold');
      
lgd = legend('show', 'Location', 'best', 'FontSize', 6, ...
             'Box', 'off', 'Orientation', 'vertical', ...
             'FontName', 'Times New Roman');
grid on;
set(gca, 'FontSize', 7, 'LineWidth', 1); 
set(gca, 'LooseInset', get(gca, 'TightInset'));
hold off;