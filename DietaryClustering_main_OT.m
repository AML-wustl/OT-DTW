% Dietary Temporal Pattern Clustering
% Entropic regularized OT-means only

clc; clear all;
%% Data Readin
data_dir = 'C:\Users\jason\Box Sync\Research\Public health\Collarboration with Yikyung\Data\';
dat = load(fullfile(data_dir,'energy_hourly_first.mat')); % dataset saves in variable "X"
X = dat.x;
fat_h = load(fullfile(data_dir,'fat_hourly_first.mat')); % dataset saves in variable "x"
carb_h = load(fullfile(data_dir,'carb_hourly_first.mat')); % dataset saves in variable "x"
prot_h = load(fullfile(data_dir,'prot_hourly_first.mat')); % dataset saves in variable "x"
Y_fat = fat_h.x; Y_carb = carb_h.x; Y_prot = prot_h.x;
Y_all = Y_fat + Y_carb + Y_prot;
Y_all_ratio = bsxfun(@rdivide, Y_all, sum(Y_all,2)); % calculated from sum of fat, carb and protein, different than Y
[num_sample, ts_len] = size(X);
MAX_K = 10;
% shift original data to avoid 0's
delta = 1e-6;
X_tmp = X + delta;
X_shift = bsxfun(@times, X_tmp, 1./(sum(X_tmp, 2)));

%% Kmeans_EROT
CH_OTD = zeros(MAX_K,1);
obj = zeros(MAX_K,1);
for K = 2:MAX_K
    [best_init_id,~] = kMeans_EROT_init(X_shift, K); % K-means++
    [member_id, centers, obj_opt, clusters] = kMeans_EROT(X_shift, K, best_init_id);
     % ground distance matrix
    [i,j] = meshgrid(1:ts_len);
    ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
    ground_d = ground_d / max(ground_d(:));
    eps = 1e-2;
    useGPU = false;
    tol = 1e-4;
    M = exp(-ground_d / eps);
    U = M .* ground_d;
    center_all = EROT_baryavg(X_shift,M,useGPU,tol);
    center_spread = zeros(K,1);
    for k = 1:K
        [tmp,~,~,~] = sinkhornTransportHotStartG(centers(k,:)',center_all',M,U);
        center_spread(k) = tmp;
    end
    CH_OTD(K) = CH_index(member_id, center_spread, obj_opt, clusters);
    obj(K) = obj_opt;
end

% [best_init_id,~] = kMeans_EROT_init(X_shift, K_opt); % K-means++
% [member_id, centers, obj_opt, clusters] = kMeans_EROT(X_shift, K_opt, best_init_id);
tic
K_opt = 7;
for iter = 1:10
    [member_id, centers, obj_opt, clusters] = kMeans_EROT(X_shift, K_opt);
end
toc


%% K-means_MallowOT
CH_OTC = zeros(MAX_K,1);
obj = zeros(MAX_K,1);
for K = 2:MAX_K
    [best_init_id,~] = kMeans_Mallows_init(X_shift, K); % K-means++
    [member_id, centers, obj_opt, clusters] = kMeans_Mallows(X_shift, K, best_init_id);
    center_all = Mallows_baryavg(X_shift);
    center_spread = zeros(K,1);
    for k = 1:K
        center_spread(k) = OT_Mallows(centers(k,:), center_all);
    end
    CH_OTC(K) = CH_index(member_id, center_spread, obj_opt, clusters);
    obj(K) = obj_opt;
end

% [best_init_id,~] = kMeans_Mallows_init(X_shift, K_opt); % K-means++
% [member_id, centers, obj_opt, clusters] = kMeans_Mallows(X_shift, K_opt, best_init_id);
tic
K_opt = 7;
for iter = 1:10
    [member_id, centers, obj_opt, clusters] = kMeans_Mallows(X_shift, K_opt);
end
toc

%% Plot J as a function of K.
h = figure(2);
width = 10;     % Width in inches
height = 50;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');

plot(2:MAX_K, CH_OTD(2:MAX_K),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
title('J as a function of K','FontSize',fsz);
xlabel('Number of cluster K','FontSize',fsz);
ylabel('Objective J','FontSize',fsz)
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\JvsK_OTD.pdf');

%% Visualization of some random samples
figure(2);
x_axis = 0:23;
width = 10;     % Width in inches
height = 75;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
sample_id = find(member_id==7);

subplot(3,2,1);
plot(x_axis, X_shift(sample_id(1),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 1','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
subplot(3,2,2);
plot(x_axis, X_shift(sample_id(2),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 2','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
subplot(3,2,3);
plot(x_axis, X_shift(sample_id(3),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 3','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
subplot(3,2,4);
plot(x_axis, X_shift(sample_id(4),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 4','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
subplot(3,2,5);
plot(x_axis, X_shift(sample_id(5),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 5','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
subplot(3,2,6);
plot(x_axis, X_shift(sample_id(6),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz);
xlim([0 24]); ylim([0 1]);
title('Sample 6','FontSize',fsz);
grid on;
set(gca, 'fontsize', fsz, 'LineWidth', alw, 'YTick',[0 0.2 0.4 0.6 0.8 1]); %<- Set properties
xlabel('Time of Eating Event','FontSize',fsz);
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size

%% Visualization of OTD-means centroids
h = figure(2);
x_axis = 0:23;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 32;      % Fontsize
lw = 5;      % LineWidth
msz = 18;       % MarkerSize
pos = get(gcf, 'Position');
p1 = plot(x_axis, centers(1,:),'-*','LineWidth',lw,'MarkerSize',msz); hold on;
p2 = plot(x_axis, centers(2,:),'-o','LineWidth',lw,'MarkerSize',msz); hold on;
p3 = plot(x_axis, centers(3,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p4 = plot(x_axis, centers(4,:),'-x','LineWidth',lw,'MarkerSize',msz); hold on;
p5 = plot(x_axis, centers(5,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p6 = plot(x_axis, centers(6,:),'-+','LineWidth',lw,'MarkerSize',msz); hold on;
p7 = plot(x_axis, centers(7,:),'-+','LineWidth',lw,'MarkerSize',msz); hold off;
% plot(x_axis, centers(7,:),'-+','Color',[0 0 0],'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(8,:),'-+','Color',[0 0 0],'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(9,:),'-+','Color',[0 0 0],'LineWidth',lw,'MarkerSize',msz); hold on;
% plot(x_axis, centers(5,:),'-+','Color',[0 0 0],'LineWidth',lw,'MarkerSize',msz); hold off;

legend([p1 p2 p3 p4 p5 p6 p7], 'C1:105','C2:287','C3:225','C4:53','C5:42','C6:156','C7:153');
xlim([0 24]); ylim([0 0.2]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily energy distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\centroids_OTD_k7.pdf');

%% Visualization of OTC-means centroid histogram
alw = 3;    % AxesLineWidth
fsz = 32;      % Fontsize
h = figure(3);
subplot(2,2,1);
p = histogram('BinEdges',0:24,'BinCounts',centers(1,:));
legend(p, 'C1:211');
xlim([0 24]); 
ylim([0 0.15]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,2);
p = histogram('BinEdges',0:24,'BinCounts',centers(2,:));
legend(p, 'C2:140');
xlim([0 24]); 
% ylim([0 0.3]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,3);
p = histogram('BinEdges',0:24,'BinCounts',centers(3,:));
legend(p, 'C3:171');
xlim([0 24]); 
ylim([0 0.15]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,4);
p = histogram('BinEdges',0:24,'BinCounts',centers(4,:));
legend(p, 'C4:234');
xlim([0 24]); 
% ylim([0 0.3]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

h = figure(4);
subplot(2,2,1);
p = histogram('BinEdges',0:24,'BinCounts',centers(5,:));
legend(p, 'C5:221');
xlim([0 24]); 
% ylim([0 0.3]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,2);
p = histogram('BinEdges',0:24,'BinCounts',centers(6,:));
legend(p, 'C6:38');
xlim([0 24]); 
% ylim([0 0.3]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,3);
p = histogram('BinEdges',0:24,'BinCounts',centers(7,:));
legend(p, 'C7:6');
xlim([0 24]); 
ylim([0 0.35]);
xlabel('Time of Eating Event','FontSize',fsz);
ylabel('Daily distribution','FontSize',fsz);
grid on;
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\centroids_OTC_k7_2.pdf');
