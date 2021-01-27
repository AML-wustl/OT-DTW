%*******************************************************************************
 % Copyright (C) 2018 Liang Wang
 % 
 % This program is free software: you can redistribute it and/or modify
 % it under the terms of the GNU General Public License as published by
 % the Free Software Foundation, version 3 of the License.
 % 
 % This program is distributed in the hope that it will be useful,
 % but WITHOUT ANY WARRANTY; without even the implied warranty of
 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 % GNU General Public License for more details.well
 % 
 % You should have received a copy of the GNU General Public License
 % along with this program.  If not, see <http://www.gnu.org/licenses/>.
 %*****************************************************************************

% OT-DTW two-stage nested Clustering project
% This main code implements our two-stage OT-DTW clustering algorithm

clc; clear all;
%%%%% Data Readin %%%%%
%% NHANES Data readin
data_dir = 'C:\Users\jason\Box Sync\Research\Public health\TemporalDietPattern\Data\NHANES\nhanesdayone\';
dat = load(fullfile(data_dir,'energy_hourly_nhanesdayone.mat')); % dataset saves in variable "X"
id = dat.x(:,1); % save the id of the person
X = dat.x(:,2:25);

[num_sample, ts_len] = size(X); % 11632 => 11631
K_outter = 7; % number of clusters
% shift original data to avoid 0's, and renormalize to get sum of 1
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2)));
% create the output matrix
output_idx = zeros(num_sample, 3);
output_idx(:,1) = id;

%% IDATA Data readin
data_dir = 'D:\Box Sync\Research\Public health\TemporalDietPattern\Data\IDATA\processed_data_old\';
dat = load(fullfile(data_dir,'energy_hourly_raw_first.mat')); % dataset saves in variable "X"
X = dat.x;
% X = str2double(dat.x(:,2:25)); % convert cell val to double val
% id = cell2mat(dat.x(:,1));
% fat_h = load(fullfile(data_dir,'fat_hourly_first.mat')); % dataset saves in variable "x"
% carb_h = load(fullfile(data_dir,'carb_hourly_first.mat')); % dataset saves in variable "x"
% prot_h = load(fullfile(data_dir,'prot_hourly_first.mat')); % dataset saves in variable "x"
% Y_fat = fat_h.x; Y_carb = carb_h.x; Y_prot = prot_h.x;
% Y_all = Y_fat + Y_carb + Y_prot;
% Y_all_ratio = bsxfun(@rdivide, Y_all, sum(Y_all,2)); % calculated from sum of fat, carb and protein, different than Y
% X = prot_h.x;
[num_sample, ts_len] = size(X); % 1021 samples
% K_outter = 4; % number of clusters
% shift original data to avoid 0's, and renormalize to get sum of 1
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2))); % the input target
% create the output matrix
output_idx = zeros(num_sample, 3);
% output_idx(:,1) = 1:1021;

%% Online Retail Data readin
data_dir = 'D:\Data\UCI ML repo\Online Retail\';
dat = csvread(fullfile(data_dir,'Sales_oneyear.csv'), 1, 1); % dataset saves in variable "X"
X = dat;
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
% output_idx(:,1) = id;

%% Smart Meter in London Data readin
data_dir = 'D:\Data\SmartMetersInLondon\';
dat = csvread(fullfile(data_dir,'all_12012013.csv'), 1, 1); % dataset saves in variable "X"
X = dat;
% fat_h = load(fullfile(data_dir,'fat_hourly_first.mat')); % dataset saves in variable "x"
% carb_h = load(fullfile(data_dir,'carb_hourly_first.mat')); % dataset saves in variable "x"
% prot_h = load(fullfile(data_dir,'prot_hourly_first.mat')); % dataset saves in variable "x"
% Y_fat = fat_h.x; Y_carb = carb_h.x; Y_prot = prot_h.x;
% Y_all = Y_fat + Y_carb + Y_prot;
% Y_all_ratio = bsxfun(@rdivide, Y_all, sum(Y_all,2)); % calculated from sum of fat, carb and protein, different than Y
[num_sample, ts_len] = size(X); 
K_outter = 7; % number of clusters
% shift original data to avoid 0's, and renormalize to get sum of 1
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2))); % the input target
% create the output matrix
output_idx = zeros(num_sample, 3);
% output_idx(:,1) = id;

%%%%% Two-stage OT-DTW algorithm %%%%%
%% 1st stage Entropy regularized OT-means
% OT parameters preparation
% ground distance matrix
[i,j] = meshgrid(1:ts_len);
ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
ground_d = ground_d / max(ground_d(:));
% barycenter parameters
useGPU = false;
eps = 1e-2;
tol = 1e-3;
M = exp(-ground_d / eps);

% choose the optimal K
minK = 2; maxK = 10;
[obj, DB, CH] = optim_K(X_shift, ground_d, useGPU, M, tol, minK, maxK, 1);

% Visualization of DB and CH index with varying K
h = figure(1);
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 3;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');
x_axis = 2:10;
p1 = plot(x_axis, CH,'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlabel('The number of cluster K','FontSize',fsz);
ylabel('OT-means CH-index','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\OnlineRetail\OT_CH.pdf');

% run the OT-means clustering with the optimal K
% K_opt_OT is determined from obj, DB and CH.
K_opt_OT = 4;
[best_init_id_EROT,~] = kMeans_EROT_init(X_shift, K_opt_OT); % K-means++
[member_id_OT, centers_OT, obj_OT, clusters_OT] = kMeans_EROT(X_shift, K_opt_OT, best_init_id_EROT);
% save the OT clustering index output
output_idx(:,2) = member_id_OT;

% Fuzzy Cmeans with OT
K_opt_OT = 4;
[best_init_id_EROT,~] = kMeans_EROT_init(X_shift, K_opt_OT); % K-means++
[~, ScoreMat_OT, centers_OT, obj_opt] = FuzzyCmeans_EROT(X_shift, K_opt_OT, best_init_id_EROT); % dim: N*K
[~, output_idx(:,2)] = max(ScoreMat_OT, [], 2);
OT_cluster_num = [sum(output_idx(:,2) == 1), sum(output_idx(:,2) == 2), sum(output_idx(:,2) == 3), sum(output_idx(:,2) == 4)];


%% 2nd stage DTW-means
OT_clusterid = 4;
tmp_idx = find(member_id_OT == OT_clusterid);
dat_to_cluster = X_shift(tmp_idx,:);

% choose the optimal K
minK = 2; maxK = 6;
[obj, DB, CH] = optim_K(dat_to_cluster, ground_d, useGPU, M, tol, minK, maxK, 2);

% run the DTW-means clustering with the optimal K
K_opt_DTW = 3;
[best_init_id_DTW,~] = kMeans_DTW_init(dat_to_cluster, K_opt_DTW); 
[member_id_DTW, centers_DTW, obj_DTW, clusters_DTW] = kMeans_DTW(dat_to_cluster, K_opt_DTW, best_init_id_DTW);

% save the DTW clustering index output
output_idx(tmp_idx,3) = member_id_DTW + 10 * OT_clusterid;

%% 2nd-stage Fuzzy DTW-means
OT_clusterid = 4;
tmp_idx = find(output_idx(:,2) == OT_clusterid);
dat_to_cluster = X_shift(tmp_idx,:);
tmp_id = id(tmp_idx,:);

% choose the optimal K
minK = 2; maxK = 6;
[obj, DB, CH] = optim_K(dat_to_cluster, ground_d, useGPU, M, tol, minK, maxK, 2);

% run the Fuzzy DTW-means clustering with the optimal K
K_opt_DTW = 4;
[best_init_id_DTW,~] = kMeans_DTW_init(dat_to_cluster, K_opt_DTW); 
[~, ScoreMat_DTW, centers_DTW, obj_opt] = FuzzyCmeans_DTW(dat_to_cluster,K_opt_DTW, best_init_id_DTW); % dim: N*K
[~, DTW_clus] = max(ScoreMat_DTW, [], 2);
output_idx(tmp_idx,3) = DTW_clus + 10 * output_idx(tmp_idx,2);
DTW_cluster_num = [sum(DTW_clus == 1), sum(DTW_clus == 2), sum(DTW_clus == 3), sum(DTW_clus == 4)];

%% Other clustering methods
% DTW_Window = length; % max value = time series length
% tic;
% [member_id, centers, clusters] = kMeans_DTW(x, K, DTW_Window);
% % [member_id, centers, clusters] = kDBA(x, K, DTW_Window);
% % [member_id, centers, clusters] = kShape(x, K);
% % [member_id, centers, clusters] = KSC(x, K);
% toc;

% Locally regularized OT + K-means
% tic;
% [member_id, centers, clusters] = kMeans_LOROT(Y, K);
% toc;

% K-medoids clustering
% tic;
% DistanceIndex = 3;    % ED=1,NCCc=2,cDTW=3, LOROT_Wasserstein=4
% DM = DMComputation(X, DistanceIndex);
% medoids = randsample(num_sample, K); % randomly pick K existing data samples as medoids.
% [member_id, new_medoids, cost, clusters] = PartitioningAroundMedoids(medoids,DM); % PAM clustering
% toc;

% save('C:\Users\jason\Box Sync\Research\Public health\Collarboration with Yikyung\Data\clusters_info_PAM_DTW.mat','member_id','new_medoids','clusters');


%% Visualization of some random samples
h = figure(1);
x_axis = 0:23;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
% sample_id = find(member_id==4);

p1 = plot(x_axis, X(1,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); hold on;
% p2 = plot(x_axis, centers_OT(2,:),'-o','Color',[0.85 0.325 0.098],'LineWidth',lw,'MarkerSize',msz); hold on;
% p3 = plot(x_axis, centers_OT(3,:),'-+','Color',[0.929 0.694 0.125],'LineWidth',lw,'MarkerSize',msz); hold on;
% p4 = plot(x_axis, centers_OT(4,:),'-x','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold on;
% legend([p1 p2 p3 p4], 'C1:53','C2:299','C3:262','C4:407');
xlim([0 24]); 
% ylim([0 0.2]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Raw Value','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

% Plot in histogram form
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
h = figure(3);
p = histogram('BinEdges',0:24,'BinCounts',X_shift(2,:),'FaceColor',[0 0.447 0.741]);
xlim([0 24]); 
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Ratio','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties


%% Visualization of 1st stage OT-means centroids
h = figure(1);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');

p1 = plot(x_axis, centers_OT(1,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers_OT(2,:),'-o','Color',[0.85 0.325 0.098],'LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers_OT(3,:),'-+','Color',[0.929 0.694 0.125],'LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers_OT(4,:),'-x','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold on;
% p5 = plot(x_axis, centers_OT(5,:),'--','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold off;
% p6 = plot(x_axis, centers_OT(6,:),'-o','LineWidth',lw,'MarkerSize',msz); hold off;
% p7 = plot(x_axis, centers_OT(7,:),'-x','LineWidth',lw,'MarkerSize',msz); hold off;

legend([p1 p2 p3 p4], 'C1:308','C2:350','C3:156','C4:207');
xlim([0 24]); 
% ylim([0 0.2]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Energy Ratio for Averaged','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\OnlineRetail\centers_OT_k6.pdf');
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\centers_OT_k4.pdf');
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\SmartMeterLondon\centers_OT_k5.pdf');
% print_to_pdf(h, 'D:\Box Sync\Research\Public Health\TemporalDietPattern\Report\images\centers_OT_weekend_k4.pdf');

%% Visualization of 2nd stage DTW-means centroids
h = figure(2);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers_DTW(1,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers_DTW(2,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers_DTW(3,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers_DTW(4,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
% p5 = plot(x_axis, centers_DTW(5,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
% p8 = plot(x_axis, centers_DTW(8,:),'-x','LineWidth',lw,'MarkerSize',msz); hold on;
% p9 = plot(x_axis, centers_DTW(9,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
% p10 = plot(x_axis, centers_DTW(10,:),'--','LineWidth',lw,'MarkerSize',msz); hold off;
legend([p1 p2 p3 p4], 'C41:101','C42:43','C43:44','C44:50');
% legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10], 'C1:9','C2:334','C3:156','C4:138','C5:180',...
%     'C6:50','C7:170','C8:131','C9:76','C10:159','location','northeastoutside');
xlim([0 24]); 
% ylim([0 1]);
xlabel('Time (hour)','FontSize',fsz);
ylabel('Energy Ratio Averaged','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\OnlineRetail\centers_DTW_1_k5.pdf');
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\centers_DTW_2_k6.pdf');
% print_to_pdf(h,'C:\Users\jason\Desktop\tkde paper\SmartMeterLondon\centers_DTW_5_k4.pdf');
% print_to_pdf(h, 'D:\Box Sync\Research\Public Health\TemporalDietPattern\Report\images\centers_DTW_weekend_2.pdf');

%% write clustering output to excel
data_dir = 'D:\Box Sync\Research\Public Health\TemporalDietPattern\Report\';
filename = 'cluster_id_042419.xlsx';
xlswrite(fullfile(data_dir,filename), id.x, 8);

%% Visualization of nutrient compositions
h = figure(1);
x_axis = 0:23;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, Y_all_ratio(1,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); hold on;
% p2 = plot(x_axis, Y_carb(8,:),'--','Color',[0.85 0.325 0.098],'LineWidth',lw,'MarkerSize',msz); hold on;
% p3 = plot(x_axis, Y_fat(8,:),'--','Color',[0.929 0.694 0.125],'LineWidth',lw,'MarkerSize',msz); hold on;
% p4 = plot(x_axis, Y_prot(8,:),'--','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold off;
xlabel('Time of Eating Event (hour)','FontSize',fsz);
ylabel('Normalized Ratio','FontSize',fsz);
% legend([p1 p2 p3 p4], 'Total intake', 'Carbohydrate intake', 'Fat intake', 'Protein intake');
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% print_to_pdf(h,'C:\Users\jason\Desktop\KDD paper\Hybrid results\intro_sample2.pdf');

%% plot 1st largest eating event of each cluster
% [~,I] = max(x,[],2);
% max_of_cluster = zeros(4,24);
% for time_id = 1:24
%     for clus_id = 1:4
%         chosen_id = find(member_id == clus_id);
%         max_of_cluster(clus_id,time_id) = sum(I(chosen_id) == time_id) / clusters(clus_id);
%     end
% end
% figure(2);
% x_axis = 0:23;
% width = 6;     % Width in inches
% height = 4.5;    % Height in inches
% alw = 0.75;    % AxesLineWidth
% fsz = 36;      % Fontsize
% lw = 5;      % LineWidth
% msz = 20;       % MarkerSize
% pos = get(gcf, 'Position');
% plot(x_axis, max_of_cluster(1,:), '-*','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max_of_cluster(2,:), '-o','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max_of_cluster(3,:), '-+','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max_of_cluster(4,:), '-x','LineWidth',lw,'MarkerSize',msz); hold on;
% hold off;
% legend('C1','C2','C3','C4');
% xlim([0 23]); ylim([0 0.35]);
% xlabel('Time of the largest Eating Event (in Hours)','FontSize',fsz);
% ylabel('Probability','FontSize',fsz);
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

%% plot 2nd largest eating event of each cluster
% [~,I_all] = sort(x,2);
% I = I_all(:,23);
% max2_of_cluster = zeros(4,24);
% for time_id = 1:24
%     for clus_id = 1:4
%         chosen_id = find(member_id == clus_id);
%         max2_of_cluster(clus_id,time_id) = sum(I(chosen_id) == time_id) / clusters(clus_id);
%     end
% end
% figure(3);
% x_axis = 0:23;
% width = 6;     % Width in inches
% height = 4.5;    % Height in inches
% alw = 0.75;    % AxesLineWidth
% fsz = 36;      % Fontsize
% lw = 5;      % LineWidth
% msz = 20;       % MarkerSize
% pos = get(gcf, 'Position');
% plot(x_axis, max2_of_cluster(1,:), '-*','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max2_of_cluster(2,:), '-o','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max2_of_cluster(3,:), '-+','LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, max2_of_cluster(4,:), '-x','LineWidth',lw,'MarkerSize',msz); hold on;
% hold off;
% legend('C1','C2','C3','C4');
% xlim([0 23]); ylim([0 0.35]);
% xlabel('Time of the 2nd largest Eating Event (in Hours)','FontSize',fsz);
% ylabel('Probability','FontSize',fsz);
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
