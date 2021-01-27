%%% Compare DBA and Wasserstein Barycenter
% test of DBA
% sequences = rand(100,20);
% mean = DBA(sequences);
% plot(mean);

%% self-made examples
% a1 = zeros(1,24); a1(8) = 0.10494341; a1(12) = 0.21046; a1(18) = 0.5; a1(18) = 0.14649;
% a1(19) = 0.1189; a1(23) = 0.4192;
% a2 = zeros(1,24); a2(8) = 0.172496293; a2(13) = 0.014316; a2(14) = 0.1822; a2(18) = 0.09418;
% a2(20) = 0.17717; a2(22) = 0.3596;
a1 = zeros(1,24); a1(8) = 0.2; a1(12) = 0.3; a1(18) = 0.5;
a2 = zeros(1,24); a2(7) = 0.1; a2(8) = 0.1; a2(11) = 0.15; a2(12) = 0.15; a2(17) = 0.25; a2(18) = 0.25;
A = [a1;a2];
% shift original data to avoid 0's, and renormalize to get sum of 1
delta = 1e-6;
A_tmp = A + delta;
A_shift = bsxfun(@times, A_tmp, 1./(sum(A_tmp, 2)));

[~,len] = size(A_shift);
% A = X_shift(3:4,:);
% cur_center = zeros(1,24);
% DTW barycenter
DBA_avg = DBA(A_shift);
% u = [0.3, 0.5, 0.1];
% DBA_avg_2 = DBA_fuzzy(A,u,2);
% OT barycenter
[i,j] = meshgrid(1:len);
ground_d = (i-j).^2; % (i-j).^2 or abs(i-j)
ground_d = ground_d / max(ground_d(:));
% barycenter parameters
eps = 1e-2;
useGPU = false;
tol = 1e-3;
M = exp(-ground_d / eps);
EROT_avg = EROT_baryavg(A_shift,M,useGPU,tol);

%% plot
h = figure(1);
x_axis = 1:24;
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 1.5;    % AxesLineWidth
fsz = 48;      % Fontsize
lw = 10;      % LineWidth
msz = 36;       % MarkerSize
pos = get(gcf, 'Position');

p1 = plot(x_axis, A(1,:),'-*','Color',[0 0.447 0.741],'LineWidth',lw,'MarkerSize',msz); grid on;
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Ratio','FontSize',fsz);
% legend('Sample 1','Sample 2');
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
% print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\barycenter_1.pdf');

h = figure(2);
plot((0:23),DBA_avg,'Color',[0 0.447 0.741],'LineWidth',lw); grid on;
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Ratio','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\barycenter_2.pdf');

h = figure(3);
plot((0:23),EROT_avg,'Color',[0 0.447 0.741],'LineWidth',lw); grid on;
xlabel('Time (hour)','FontSize',fsz);
ylabel('Normalized Ratio','FontSize',fsz);
set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties
print_to_pdf(h,'C:\Users\jason\Box Sync\Research\My paper\DietaryTimeSeriesClustering\IEEE TKDE\plots\barycenter_3.pdf');

% h = figure(4);
% plot((0:23),DBA_avg,'Color',[0 0.447 0.741],'LineWidth',lw); hold on;
% plot((0:23),DBA_avg_2,'Color',[0.85 0.325 0.098],'LineWidth',lw); hold off;
% grid on;
% xlabel('Time of Eating Event (hour)','FontSize',fsz);
% ylabel('Normalized Ratio','FontSize',fsz);
% legend('DTW barycenter');
% set(gca, 'fontsize', fsz, 'LineWidth', alw); %<- Set properties