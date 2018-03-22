function [solutionTime, displ_mean, displ_mean_rate, vel_mean, vel_mean_rate, strE_tot, strE_tot_rate] = plotConvergence(tarDIR,bodyFN,param)

% Load data
matDispl = fullfile(tarDIR,strcat(bodyFN,'_LoadStepDisplVel.mat'));
matStrE = fullfile(tarDIR,strcat(bodyFN,'_LoadStepStrE.mat'));
load(matDispl,'displ_vel','solutionTime');
load(matStrE,'strE','trussStrE','misalignHJE','misalignDSDNA','solutionTime_outfile');

% Check data
assert(max(abs((solutionTime(2:end)-solutionTime_outfile)./solutionTime_outfile)) < 1e-4);
assert(size(displ_vel,2)==6);
nLoadStep = numel(solutionTime);
assert(size(displ_vel,3)==nLoadStep && size(strE,2)==nLoadStep);

% Time history of the mean nodal displacement
displ_mean = zeros(nLoadStep,1);
for i = 1:nLoadStep
    displ_mean(i) = mean(sqrt(sum(displ_vel(:,1:3,i).^2,2)));
end
displ_mean_rate = timeRate(solutionTime,displ_mean);

% Time history of the mean nodal velocity
vel_mean = zeros(nLoadStep,1);
for i = 1:nLoadStep
    vel_mean(i) = mean(sqrt(sum(displ_vel(:,4:6,i).^2,2)));
end
vel_mean_rate = timeRate(solutionTime,vel_mean);

% Time history of the total strain energy
strE_tot = zeros(nLoadStep,1);
for i = 1:nLoadStep
    if(i==1)
        strE_tot(i) = sum(strE(:,i)) + sum(trussStrE(:,i));
    else
        strE_tot(i) = sum(strE(:,i)) + sum(trussStrE(:,i)) + sum(alignE(misalignHJE(:,:,i-1),param.alignHJE))/param.KbT + sum(alignE(misalignDSDNA(:,:,i-1),param.alignDSDNA))/param.KbT;
    end
end
strE_tot_rate = timeRate(solutionTime,strE_tot);


matStrE = fullfile(tarDIR,strcat(bodyFN,'_Conv.mat'));
save(matStrE,'strE_tot');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
conv_crit = 1e-3;
stepmax = 110;
xmin = 1;
xmax1 = 1e10;
xmax2 = 200;
ymin = 1e-6;
ymax = 1e5;
figWidth = 6.5;
figHeight = 3;
pos1x = 0.08;
pos2x = 0.5+pos1x;
pos1y = 0.15;
pos2y = 0.5+pos1y;
pos3 = 0.5-pos1x-0.02;
pos4 = 0.5-pos1y-0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(101);
set(gcf, 'Units','inches','PaperUnits','inches');
pos = get(gcf, 'Position');
pos(3) = figWidth;
pos(4) = figHeight;
set(gcf, 'Position',pos);
set(gcf,'PaperSize',[figWidth figHeight], 'PaperPosition',[0 0 figWidth figHeight]);

subplot(221);
semilogx(solutionTime,displ_mean,'.'); hold on;
semilogx(solutionTime,displ_mean,'-');
set(gca,'XLim',[xmin xmax1], 'FontSize',6);
xlabel('time','FontSize',6); ylabel(strcat('mean nodal displ. mag. (', char(197), ')'), 'FontSize',6);
set(gca,'Position',[pos1x pos2y pos3 pos4]);

subplot(222);
plot(solutionTime(1:stepmax),displ_mean(1:stepmax),'.'); hold on;
plot(solutionTime(1:stepmax),displ_mean(1:stepmax),'-');
set(gca,'XLim',[xmin xmax2], 'FontSize',6);
xlabel('time','FontSize',6); ylabel(strcat('mean nodal displ. mag. (', char(197), ')'), 'FontSize',6);
set(gca,'Position',[pos2x pos2y pos3 pos4]);

subplot(223);
loglog(solutionTime,vel_mean,'.'); hold on;
loglog(solutionTime,vel_mean,'-'); hold on;
loglog([xmin xmax1],[conv_crit conv_crit], ':k');
set(gca,'XLim',[xmin xmax1], 'YLim',[ymin ymax], 'FontSize',6);
xlabel('time','FontSize',6); ylabel(strcat('mean nodal vel. mag. (', char(197), '/unit time)'), 'FontSize',6);
set(gca,'Position',[pos1x pos1y pos3 pos4]);

subplot(224);
semilogy(solutionTime(1:stepmax),vel_mean(1:stepmax),'.'); hold on;
semilogy(solutionTime(1:stepmax),vel_mean(1:stepmax),'-'); hold on;
semilogy([xmin xmax2],[conv_crit conv_crit], ':k');
set(gca,'XLim',[xmin xmax2], 'YLim',[ymin ymax], 'FontSize',6);
xlabel('time','FontSize',6); ylabel(strcat('mean nodal vel. mag. (', char(197), '/unit time)'), 'FontSize',6);
set(gca,'Position',[pos2x pos1y pos3 pos4]);

print(gcf, '-dtiff', '-r300', fullfile(tarDIR,strcat(bodyFN,'_LoadStepDispl.tif')));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: strain energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(201);
set(gcf, 'Units','inches','PaperUnits','inches');
pos = get(gcf, 'Position');
pos(3) = figWidth;
pos(4) = figHeight;
set(gcf, 'Position',pos);
set(gcf,'PaperSize',[figWidth figHeight], 'PaperPosition',[0 0 figWidth figHeight]);

subplot(221);
semilogx(solutionTime,strE_tot,'.'); hold on;
semilogx(solutionTime,strE_tot,'-');
set(gca,'XLim',[xmin xmax1], 'FontSize',6);
xlabel('time','FontSize',6); ylabel('total strain energy (k_BT)', 'FontSize',6);
set(gca,'Position',[pos1x pos2y pos3 pos4]);

subplot(222);
plot(solutionTime(1:stepmax),strE_tot(1:stepmax),'.'); hold on;
plot(solutionTime(1:stepmax),strE_tot(1:stepmax),'-');
set(gca,'XLim',[xmin xmax2], 'FontSize',6);
xlabel('time','FontSize',6); ylabel('total strain energy (k_BT)', 'FontSize',6);
set(gca,'Position',[pos2x pos2y pos3 pos4]);

subplot(223);
loglog(solutionTime,strE_tot_rate,'.'); hold on;
loglog(solutionTime,strE_tot_rate,'-'); hold on;
set(gca,'XLim',[xmin xmax1], 'YLim',[ymin ymax], 'FontSize',6);
xlabel('time','FontSize',6); ylabel('rate of total strain energy (k_BT/unit time)', 'FontSize',6);
set(gca,'Position',[pos1x pos1y pos3 pos4]);

subplot(224);
semilogy(solutionTime(1:stepmax),strE_tot_rate(1:stepmax),'.'); hold on;
semilogy(solutionTime(1:stepmax),strE_tot_rate(1:stepmax),'-'); hold on;
set(gca,'XLim',[xmin xmax2], 'YLim',[ymin ymax], 'FontSize',6);
xlabel('time','FontSize',6); ylabel('rate of total strain energy (k_BT/unit time)', 'FontSize',6);
set(gca,'Position',[pos2x pos1y pos3 pos4]);

print(gcf, '-dtiff', '-r300', fullfile(tarDIR,strcat(bodyFN,'_LoadStepStrE.tif')));

end


% Mechanical energy for alignment elements
function [E,E_component] = alignE(misalign,K)

dof = zeros(6,size(misalign,2));
dof(1:3,:) = misalign(2:4,:);
dof(4:6,:) = bsxfun(@times,misalign(5,:),misalign(6:8,:));
E_component = 1/2 * (K.KTC1*dof(1,:).^2 + K.KTC2*dof(2,:).^2 + K.KTC3*dof(3,:).^2 + ...
                     K.KRC1*dof(4,:).^2 + K.KRC2*dof(5,:).^2 + K.KRC3*dof(6,:).^2);
E = sum(E_component,1);

end


% Calculate time rate
function r = timeRate(t,x)

assert(iscolumn(t) && iscolumn(x) && numel(t)==numel(x));
r = zeros(size(t));
r(1) = abs((x(2)-x(1))/(t(2)-t(1)));
r(end) = abs((x(end)-x(end-1))/(t(end)-t(end-1)));

x_prev = x(1:end-2);
x_curr = x(2:end-1);
x_next = x(3:end);
t_prev = t(1:end-2);
t_curr = t(2:end-1);
t_next = t(3:end);

r_prev = (x_curr-x_prev)./(t_curr-t_prev);
r_next = (x_next-x_curr)./(t_next-t_curr);

r(2:end-1) = abs(mean([r_prev r_next],2));

end