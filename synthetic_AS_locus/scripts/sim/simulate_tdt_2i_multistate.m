clear all; clc;
%% Load data

opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Time", "GFP", "Tomato", "VarName4"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Time", "GFP", "Tomato", "VarName4"], "DecimalSeparator", ",");

% Import the data
%2i
esc_ini1 = readtable("data_fig6_2i_ini1.txt", opts);
norm_fac_esc = mean(esc_ini1.GFP(1:5));
norm_fac_escT = mean(esc_ini1.Tomato(1:5));

esc_ini2 = readtable("data_fig6_2i_ini2.txt", opts);
norm_fac_esc2 = mean(esc_ini2.GFP(1:5));
norm_fac_esc2T = mean(esc_ini2.Tomato(1:5));

% Diff
diff_ini1 = readtable("data_fig6_diff_ini1.txt", opts);

norm_fac = mean(diff_ini1.GFP(1:3));
norm_facT = mean(diff_ini1.Tomato(1:3));

% Ini2
diff_ini2 = readtable("data_fig6_diff_ini2.txt", opts);

norm_fac2 = mean(diff_ini2.GFP(1:3));
norm_fac2T = mean(diff_ini2.Tomato(1:3));

%control data diff tdt
opts = setvaropts(opts, ["Time", "GFP", "Tomato", "VarName4"], "DecimalSeparator", ".");
diff_cntrl_t = readtable("data_diff_fullDox.csv", opts);
norm_fac2T_c = mean(diff_ini2.Tomato(1:3));

%% Fitting tdt 2i using profile likelihood

data_tdt = sortrows([esc_ini2.Time (esc_ini2.Tomato-norm_fac_escT)/(norm_fac_esc2T-norm_fac_escT)]);
[g, h] = findgroups(data_tdt(:,1)); 
data_mean_tdt = [h splitapply(@(x)mean(x,1), data_tdt(:, 2), findgroups(data_tdt(:,1)))];
data_std_tdt_all = [h splitapply(@(x)std(x), data_tdt(:, 2), findgroups(data_tdt(:,1)))];

% leave out last tp?
data_std_tdt = data_std_tdt_all(1:end,:);
t_data = data_mean_tdt(1:end,1);
protein_data = data_mean_tdt(1:end,2);


%% Settings and parameters
kdeg_tdt = log(2)/(26*0.86);
kdil = log(2)/12; % dilution due to cell division (cell cycle length ~12h (Waisman2018))
kdd_tdt = kdeg_tdt + kdil; % total rate of tdt removal due to degradation & dilution

n_range = 1:11;                      % Number of ON states to test (= intermediate states+1)
k_vals = logspace(-1.1, 0.2, 500);        % Transition rate grid (1/h)
k_deg = kdd_tdt;                     % Known degradation rate
t = 0:100;                         % Time points (h)
y0_template = @(n) [1; zeros(n, 1);1]; % Initial state: ON_1 - ON_n, OFF, protein

%% Loop over different promoter structures
results = struct();

for ni = 1:length(n_range)
    n = n_range(ni);
    y0 = y0_template(n);
    chi2_profile = zeros(size(k_vals));

    for ki = 1:length(k_vals)
        k = k_vals(ki);
        [~, y] = ode45(@(t,y) promoter_model(t, y, k, k_deg, n), t_data, y0);
        protein_sim = y(:, end);
        residuals = (protein_sim - protein_data)./data_std_tdt(:,2);
        chi2_profile(ki) = sum(residuals.^2);
    end

    % Store results
    [chi2_min, best_idx] = min(chi2_profile);
    k_best = k_vals(best_idx);
    

    results(ni).n = n;
    results(ni).k_vals = k_vals;
    results(ni).chi2 = chi2_profile;
    results(ni).k_best = k_best;
    results(ni).chi2_min = chi2_min;
    
end

%% Plot profile likelihoods
figure(1);clf; hold on;
for ni = 1:length(n_range)
    plot(log10(results(ni).k_vals), results(ni).chi2 - results(ni).chi2_min, ...
        'DisplayName', ['n = ' num2str(results(ni).n)]);
end
yline(3.84, 'r--', '95% CI');
xlabel('log_{10}(k)');
ylabel('\Delta\chi^2');
legend('show');
set(gca, 'YLim', [0,20]);
title('Profile Likelihoods for k across promoter models');
grid on;

%% Model comparison summary
fprintf('\nModel Comparison Summary:\n');
fprintf(' n   k_best     chi2\n');
summary = [];
for ni = 1:length(results)
    r = results(ni);
    fprintf('%2d   %.4f   %.2f\n', ...
        r.n, r.k_best, r.chi2_min);
    summary = [summary; r.n, r.k_best, r.chi2_min];
end
[val idx] = min(summary(:,3));
best_n = min(summary(idx,1));
best_k = min(summary(idx,2));
%% Extract 95% CI for best n
% Normalize
chi2_rel = results(idx).chi2 - min(results(idx).chi2); % Relative to best fit
[val best_idx] = min(results(idx).chi2); % extract index of best k

% Interpolate to find where chi2 crosses the CI threshold
alpha = 0.95;  % confidence level
threshold = chi2inv(alpha, 1);

% Search left of best fit
left_idx = find(chi2_rel(1:best_idx) < threshold, 1, 'first');
if ~isempty(left_idx) && left_idx < best_idx
    k_low = interp1( ...
        chi2_rel(left_idx-1:left_idx), ...
        k_vals(left_idx-1:left_idx), ...
        threshold);
else
    k_low = NaN;  % unbounded
end

% Search right of best fit
right_idx = find(chi2_rel(best_idx:end) > threshold, 1, 'first') + best_idx - 1;
if ~isempty(right_idx) && right_idx > best_idx
    k_high = interp1( ...
        chi2_rel(right_idx-1:right_idx), ...
        k_vals(right_idx-1:right_idx), ...
        threshold);
else
    k_high = NaN;  % unbounded
end

% Print results
fprintf('\n95%% Confidence Interval for k_on:\n');
fprintf('Best fit: %.4f\n', best_k);
if isnan(k_low)
    fprintf('Lower bound: unbounded\n');
else
    fprintf('Lower bound: %.4f\n', k_low);
end
if isnan(k_high)
    fprintf('Upper bound: unbounded\n');
else
    fprintf('Upper bound: %.4f\n', k_high);
end


%% Plot simulation 
tspan=[1:100];
n=best_n; k=best_k; y0=y0_template(n);
[~, y] = ode45(@(t,y) promoter_model(t, y, k, k_deg, n), tspan, y0);
% simulate wo delay
n=1; k = 1000;y0=y0_template(n);
[~, y2] = ode45(@(t,y) promoter_model(t, y, k, k_deg, n), tspan, y0);
% simulate ini1
y0_template_ini1 = @(n) [zeros(n, 1);1; 0]; % Initial state: ON_1 - ON_n, OFF, protein
n=best_n; k=best_k; y0=y0_template_ini1(n);
[~, y3] = ode45(@(t,y) promoter_model(t, y, k, k_deg, n), tspan, y0);
figure(2); clf;
hold off;
set(gcf,'DefaultAxesColorOrder',[1 0 0;0 1 0])

plot(esc_ini2.Time,log10(esc_ini2.Tomato),'o')
hold on
plot(tspan,log10(y(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT))
hold on
plot(tspan,log10(y2(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT),':', 'color','green')
hold on
plot(tspan,log10(y3(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT),':', 'color','black')

legend('data','sim w/ delay', 'sim w/o delay')
title ('ESC - Ini2 - Dox - Tomato')
xlabel('Time (h)')
ylabel ('norm. MFI (log10)')
%saveas(gcf,'Fit_ESC_Ini2_Tomato.pdf')

%%
% write simulation into csv file
writematrix([tspan',y(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT],'sim_best_fit_2i_Ini2_tdt.csv');
writematrix([tspan',y2(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT],'sim_no_delay_2i_Ini2_tdt.csv');
writematrix([tspan',y3(:,end).*(norm_fac_esc2T-norm_fac_escT)+norm_fac_escT],'sim_2i_Ini1_tdt.csv');


%% ODE function: chain of promoter states + protein

function dydt = promoter_model(t, y, k, k_deg, n)
    dydt = zeros(n+2, 1);  % n OFF states + ON + protein
    ON = y(1:n);        % promoter states from ON_1 to ON_n 
    OFF = y(n+1); % final promoter ON state
    P = y(end);            % protein

    % promoter transitions
    dydt(1) = -k * ON(1);
    for i = 2:n
        dydt(i) = k * (ON(i-1) - ON(i));
    end
    dydt(n+1) = k * ON(n);
    % Protein 
    dydt(end) = k_deg * (1-OFF - P);
end



