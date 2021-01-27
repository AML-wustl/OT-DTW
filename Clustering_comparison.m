clc; clear all;
%%%%% Data Readin %%%%%
%% IDATA Data readin
data_dir = 'C:\Users\jason\Box Sync\Research\Public health\TemporalDietPattern\Data\IDATA\processed_data\';
dat = load(fullfile(data_dir,'energy_hourly_raw_first.mat')); % dataset saves in variable "X"
X = dat.x;
% fat_h = load(fullfile(data_dir,'fat_hourly_first.mat')); % dataset saves in variable "x"
% carb_h = load(fullfile(data_dir,'carb_hourly_first.mat')); % dataset saves in variable "x"
% prot_h = load(fullfile(data_dir,'prot_hourly_first.mat')); % dataset saves in variable "x"
% Y_fat = fat_h.x; Y_carb = carb_h.x; Y_prot = prot_h.x;
% Y_all = Y_fat + Y_carb + Y_prot;
% Y_all_ratio = bsxfun(@rdivide, Y_all, sum(Y_all,2)); % calculated from sum of fat, carb and protein, different than Y
[num_sample, ts_len] = size(X); % 1021 samples
K_outter = 7; % number of clusters
% shift original data to avoid 0's, and renormalize to get sum of 1
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2))); % the input target
% create the output matrix
output_idx = zeros(num_sample, 3);

%% OT-OT
% 1st stage, run the OT-means clustering with the optimal K
% K_opt_OT is determined from obj, DB and CH.
K_opt_OT = 4;
[best_init_id_EROT,~] = kMeans_EROT_init(X_shift, K_opt_OT); % K-means++
[member_id_OT, centers_OT, obj_OT, clusters_OT] = kMeans_EROT(X_shift, K_opt_OT, best_init_id_EROT);
% save the OT clustering index output
output_idx(:,2) = member_id_OT;

% 2nd stage
OT_clusterid = 2;
tmp_idx = find(member_id_OT == OT_clusterid);
dat_to_cluster = X_shift(tmp_idx,:);

% choose the optimal K
% ground distance matrix
[i,j] = meshgrid(1:ts_len);
ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
ground_d = ground_d / max(ground_d(:));
% barycenter parameters
useGPU = false;
eps = 1e-2;
tol = 1e-3;
M = exp(-ground_d / eps);
minK = 2; maxK = 10;
[obj, DB, CH] = optim_K(dat_to_cluster, ground_d, useGPU, M, tol, minK, maxK, 1);
K_opt_2nd = 7;
[best_init_id_2nd,~] = kMeans_EROT_init(dat_to_cluster, K_opt_2nd); % K-means++
[member_id_2nd, centers_2nd, obj_2nd, clusters_2nd] = kMeans_EROT(dat_to_cluster, K_opt_2nd, best_init_id_2nd);
% save the OT clustering index output
output_idx(tmp_idx,3) = member_id_2nd + 10 * OT_clusterid;

% Visualization of 2nd stage OT-means centroids
h = figure(2);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers_2nd(1,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers_2nd(2,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers_2nd(3,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers_2nd(4,:),'-x','LineWidth',lw,'MarkerSize',msz); hold on;
p5 = plot(x_axis, centers_2nd(5,:),'--','LineWidth',lw,'MarkerSize',msz); hold on;
p6 = plot(x_axis, centers_2nd(6,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p7 = plot(x_axis, centers_2nd(7,:),'-+','LineWidth',lw,'MarkerSize',msz); hold off;
legend([p1 p2 p3 p4 p5 p6 p7], 'C1:46','C2:67','C3:42','C4:133','C5:17','C6:8','C7:6','location','northeastoutside');
xlim([0 24]); 
% ylim([0 1]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Energy Ratio','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\TemporalDietaryData\centers_OT_2nd_2_k7.pdf');

%% DTW-DTW
% 1st stage, run the OT-means clustering with the optimal K
% K_opt_OT is determined from obj, DB and CH.
% choose the optimal K
minK = 2; maxK = 10;
[obj, DB, CH] = optim_K(X_shift, ground_d, useGPU, M, tol, minK, maxK, 2);
K_opt_1st = 4;
[best_init_id_1st,~] = kMeans_DTW_init(X_shift, K_opt_1st); % K-means++
[member_id_1st, centers_1st, obj_1st, clusters_1st] = kMeans_DTW(X_shift, K_opt_1st, best_init_id_1st);
% save the OT clustering index output
output_idx(:,2) = member_id_1st;

% 2nd stage
clusterid = 4;
tmp_idx = find(member_id_1st == clusterid);
dat_to_cluster = X_shift(tmp_idx,:);
% choose the optimal K
minK = 2; maxK = 10;
[obj, DB, CH] = optim_K(dat_to_cluster, ground_d, useGPU, M, tol, minK, maxK, 2);

K_opt_2nd = 5;
[best_init_id_2nd,~] = kMeans_DTW_init(dat_to_cluster, K_opt_2nd); % K-means++
[member_id_2nd, centers_2nd, obj_2nd, clusters_2nd] = kMeans_DTW(dat_to_cluster, K_opt_2nd, best_init_id_2nd);
% save the OT clustering index output
output_idx(tmp_idx,3) = member_id_2nd + 10 * OT_clusterid;

% Visualization of centroids
h = figure(2);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers_2nd(1,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers_2nd(2,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers_2nd(3,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers_2nd(4,:),'-x','LineWidth',lw,'MarkerSize',msz); hold on;
p5 = plot(x_axis, centers_2nd(5,:),'--','LineWidth',lw,'MarkerSize',msz); hold on;
% p6 = plot(x_axis, centers_2nd(6,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
% p7 = plot(x_axis, centers_2nd(7,:),'-+','LineWidth',lw,'MarkerSize',msz); hold off;

legend([p1 p2 p3 p4 p5], 'C1:44','C2:83','C3:66','C4:98','C5:6','location','northeastoutside');
xlim([0 24]); 
% ylim([0 1]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Energy Ratio','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\TemporalDietaryData\centers_DTW_2nd_4_k5.pdf');

%% MDTW-MDTW
% 1st stage, run the MDTW-means clustering with the optimal K
% K_opt_1st is determined from obj, DB and CH.
K_opt_1st = 4;
% [best_init_id_1st,~] = kMeans_DTW_init(X_shift, K_opt_1st); % K-means++
[member_id_1st, centers_1st, obj_1st, clusters_1st] = kMeans_MDTW(X_shift, K_opt_1st);
% save the OT clustering index output
output_idx(:,2) = member_id_1st;

% 2nd stage
clusterid = 4;
tmp_idx = find(member_id_1st == clusterid);
dat_to_cluster = X_shift(tmp_idx,:);
K_opt_2nd = 5;
[best_init_id_2nd,~] = kMeans_DTW_init(dat_to_cluster, K_opt_2nd); % K-means++
[member_id_2nd, centers_2nd, obj_2nd, clusters_2nd] = kMeans_DTW(dat_to_cluster, K_opt_2nd, best_init_id_2nd);
% save the OT clustering index output
output_idx(tmp_idx,3) = member_id_2nd + 10 * OT_clusterid;

% Visualization of centroids
h = figure(2);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers_1st(1,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers_1st(2,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers_1st(3,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers_1st(4,:),'-x','LineWidth',lw,'MarkerSize',msz); hold on;
% p5 = plot(x_axis, centers_2nd(5,:),'--','LineWidth',lw,'MarkerSize',msz); hold on;
% p6 = plot(x_axis, centers_2nd(6,:),'-*','LineWidth',lw,'MarkerSize',msz); hold off;

legend([p1 p2 p3 p4], 'C1:163','C2:151','C3:513','C4:194');
xlim([0 24]); 
% ylim([0 1]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Energy Consumption Ratio','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\TemporalDietaryData\centers_DTW_2nd_4_k5.pdf');