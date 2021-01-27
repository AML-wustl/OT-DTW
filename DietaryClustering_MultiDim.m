% Dietary Pattern Clustering project

clc; clear all;
%% Data Readin
data_dir = 'C:\Users\jason\Box Sync\Research\Public health\Collarboration with Yikyung\Data\';
fat_h = load(fullfile(data_dir,'fat_hourly_1.mat')); % dataset saves in variable "x"
carb_h = load(fullfile(data_dir,'carb_hourly_1.mat')); % dataset saves in variable "x"
prot_h = load(fullfile(data_dir,'prot_hourly_1.mat')); % dataset saves in variable "x"
fat_dat = fat_h.x; carb_dat = carb_h.x; prot_dat = prot_h.x;

[num_sample, length] = size(fat_dat);
K = 4; % number of clusters

%% Calculate Distance Matrix
tic;
DistanceIndex = 3;    % ED = 1,NCCc = 2,cDTW = 3
DM_f = DMComputation(fat_dat, DistanceIndex);
DM_c = DMComputation(carb_dat, DistanceIndex);
DM_p = DMComputation(prot_dat, DistanceIndex);
DM = DM_f + DM_c + DM_p;
medoids = randsample(num_sample, K);
[member_id, new_medoids, cost, clusters] = PartitioningAroundMedoids(medoids,DM_c);
toc;

%% Kmeans+DTW clustering
DTW_Window = 24; % max value = time series length
tic;
% [member_id, centers, clusters] = kMeans_ED(x, K);
% [member_id, centers, clusters] = kDBA(x, K, DTW_Window);
% [member_id, centers, clusters] = kShape(x, K);
% [member_id, centers, clusters] = KSC(x, K);
toc;
save('C:\Users\jason\Box Sync\Research\Public health\Collarboration with Yikyung\Data\clusters_multi_PAM_DTW.mat','member_id','new_medoids','clusters');

%% Plot
load(fullfile(data_dir,'energy_hourly_1_new.mat')); % dataset saves in variable "x"
load(fullfile(data_dir,'clusters_multi_PAM_DTW.mat'));
% plot centroid energy time series for visualization
figure(1);
x_axis = 0:23;
width = 6;     % Width in inches
height = 2.25;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 28;      % Fontsize
lw = 5;      % LineWidth
msz = 20;       % MarkerSize
pos = get(gcf, 'Position');
plot(x_axis, prot_dat(new_medoids(1),:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); hold on;
plot(x_axis, prot_dat(new_medoids(2),:),'-o','Color',[0.85 0.325 0.098],'LineWidth',lw,'MarkerSize',msz); hold on;
plot(x_axis, prot_dat(new_medoids(3),:),'-+','Color',[0.929 0.694 0.125],'LineWidth',lw,'MarkerSize',msz); hold on;
plot(x_axis, prot_dat(new_medoids(4),:),'-x','Color',[0.494 0.184 0.556],'LineWidth',lw,'MarkerSize',msz); hold off;
legend('C1','C2','C3','C4');
xlim([0 23]); ylim([0 0.8]);
xlabel('Time of Eating Event for centroid of 4 clusters','FontSize',fsz);
ylabel('Normalized daily distribution for protein','FontSize',fsz);
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties

% plot 1st and 2nd largest eating event of each cluster
