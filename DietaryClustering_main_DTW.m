% Dietary Temporal Pattern Clustering
% DTW-means clustering only

clc; clear all;
%% Data Readin
data_dir = 'C:\Users\jason\Box Sync\Research\Public health\Collarboration with Yikyung\Data\';
dat = load(fullfile(data_dir,'energy_hourly_1.mat')); % dataset saves in variable "X"
X = dat.x;
dat = load(fullfile(data_dir,'energy_hourly_all.mat')); % dataset saves in variable "Y"
fat_h = load(fullfile(data_dir,'fat_hourly_all.mat')); % dataset saves in variable "x"
carb_h = load(fullfile(data_dir,'carb_hourly_all.mat')); % dataset saves in variable "x"
prot_h = load(fullfile(data_dir,'prot_hourly_all.mat')); % dataset saves in variable "x"
Y_fat = fat_h.x'; Y_carb = carb_h.x'; Y_prot = prot_h.x';
Y = dat.x';
Y_all = Y_fat + Y_carb + Y_prot;
Y_all_ratio = bsxfun(@rdivide, Y_all, sum(Y_all,2)); % calculated from sum of fat, carb and protein, different than Y
[num_sample, ts_len] = size(X);
K = 18; % number of clusters
% shift original data to avoid 0's
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2)));

%% DTW + K-means
tic;
run_num = 100; % number of randomized initializations
init_id_DTW = zeros(num_sample,run_num);
member_id_DTW = zeros(num_sample,run_num);
centers_DTW = zeros(K,ts_len,run_num);
obj_DTW = zeros(run_num,1);
for i = 1:run_num
    disp(i);
    [init_id, member_id, centers, obj_opt, ~] = kMeans_DTW(X_shift, K);
    init_id_DTW(:,i) = init_id;
    member_id_DTW(:,i) = member_id;
    centers_DTW(:,:,i) = centers;
    obj_DTW(i) = obj_opt;
end
toc;
% post-processing of clustering results
% load('DTW_results.mat');
[sorted_obj_DTW,I] = sort(obj_DTW);
plot(sorted_obj_DTW);
best_init_id_DTW = init_id_DTW(:, I(1));
[init_id, member_id, centers, obj_opt, clusters] = kMeans_DTW(X_shift, K, best_init_id_DTW);
% save('DTW_results.mat','init_id_DTW','member_id_DTW','centers_DTW','obj_DTW');

%% Visualization of some random samples
h = figure(1);
x_axis = 0:23;
width = 10;     % Width in inches
height = 50;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 24;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
samples = X_shift(member_id==7,:);

subplot(5,1,1);
plot(x_axis, samples(1,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);grid on;
title('Sample 1','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
subplot(5,1,2);
plot(x_axis, samples(2,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);grid on;
title('Sample 2','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
subplot(5,1,3);
plot(x_axis, samples(3,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);grid on;
title('Sample 3','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
subplot(5,1,4);
plot(x_axis, samples(4,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);grid on;
title('Sample 4','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
subplot(5,1,5);
plot(x_axis, samples(5,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);grid on;
title('Sample 5','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
xlabel('Time of Eating Event','FontSize',fsz); grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% print_to_pdf(h,'C:\Users\jason\Desktop\KDD paper\Hybrid results\DTWonly_samples_2.pdf');

%% Visualization of K-means centroids
h = figure(2);
x_axis = 0:23;
width = 12.5;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers(12,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers(17,:),'-o','Color',[0.85 0.325 0.098],'LineWidth',lw,'MarkerSize',msz); hold on;
% p3 = plot(x_axis, centers(18,:),'-+','Color',[0.929 0.694 0.125],'LineWidth',lw,'MarkerSize',msz); hold on;
% p4 = plot(x_axis, centers(14,:),'-x','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold on;
% p5 = plot(x_axis, centers(15,:),'-x','Color',[0 1 0],'LineWidth',lw,'MarkerSize',msz); hold off;
% plot(x_axis, centers(15,:),'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(6,:),'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(7,:),'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(8,:),'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(9,:),'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(10,:),'LineWidth',lw,'MarkerSize',msz); hold off;

legend([p1 p2],'C1','C2');
xlim([0 24]); ylim([0 1]);
xlabel('Time of Eating Event for centers','FontSize',fsz);
ylabel('Normalized daily distribution','FontSize',fsz);
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
