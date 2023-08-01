% =========================================================================
% All the results
% =========================================================================
clear all; close all;
% First analysis
t1 = readtable('mr_results_scaled.csv');
% Wigthman
t2 = readtable('mr_results_2_scaled.csv');

blue = [114 172 185]./255;
or   = [239 186 0]./255;

f = figure(1); hold on; grid on; box on;
errorbar(flip(t1.OR), [1.05, 1.45], flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), [0.95, 1.35], flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or,"Color",or)
xline(1,'LineStyle','--','LineWidth',2)
ylim([.8,1.6])
ax = gca; ax.YTick = [1, 1.4]; ax.YTickLabel = {'SBP','DBP'}; ax.XLim = [0.4 5.5]; ax.XTick = [0.5,1,2,3,4,5]
l = legend('Lambert et al. 2013', 'Wightman et al. 2021','Location','Northoutside','NumColumns',2);
set(gca, 'XScale', 'log')
xlabel('OR per 100mg daily use of sildenafil')
print(f,'./Figures/All.png','-r300','-dpng')

%% ========================================================================
% Two step results // Diastolic blood pressure
% =========================================================================
clear all; close all;

% First analysis
t1 = readtable('mr_results_scaled.csv');
t1 = t1(1,:);

t2 = readtable('mr_results_twoStep_DBP_scaled.csv');
t2 = t2(1:14,:);

blue = [114 172 185]./255;

f = figure(1); hold on; grid on; box on;
xline(1,'LineStyle','--','LineWidth',2,'Color','k')
errorbar(flip(t1.OR), 15, flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), 1:14, flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
ax = gca; ax.YTick = [1:15]; ax.YTickLabel = flip({'Main','BMI','Impedance of leg (right)','Impedance of leg (left)','Impedance of arm (right)','Impedance of arm (left)','Impedance of whole body','Standing height','Plateletcrit','Myeloid white cell count','Platelet count','White blood cell count','Coronary artery disease','Granulocyte count','Sum basophil neutrophil counts'});
ax.XTick = [.5 .75 1 1.5 2];
set(gca,'XScale','log')
xlim([.45,2.3])
ylim([0.5,15.5])
fill([.45 2.3 2.3 .45],[14.5 14.5 15.5 15.5],[.7 .7 .7],'FaceAlpha',0.2,'EdgeColor',[.7 .7 .7],'EdgeAlpha',0.2)
title('Diastolic blood pressure')
xlabel('OR per 100mg daily use of sildenafil')
print(f,'./Figures/T2MR_DBP.png','-r300','-dpng')

%% ========================================================================
% Systolic blood pressure
% =========================================================================
clear all; close all;

% First analysis
t1 = readtable('mr_results_scaled.csv');
t1 = t1(2,:);

t2 = readtable('mr_results_twoStep_SBP_scaled.csv');
t2 = t2(15:end,:);

blue = [114 172 185]./255;

f = figure(1); hold on; grid on; box on;
xline(1,'LineStyle','--','LineWidth',2,'Color','k')
errorbar(flip(t1.OR), 15, flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), 1:14, flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
ax = gca; ax.YTick = [1:15]; ax.YTickLabel = flip({'Main','BMI','Impedance of leg (right)','Impedance of leg (left)','Impedance of arm (right)','Impedance of arm (left)','Impedance of whole body','Standing height','Plateletcrit','Myeloid white cell count','Platelet count','White blood cell count','Coronary artery disease','Granulocyte count','Sum basophil neutrophil counts'});
ax.XTick = [.5 1 2 3 4 5];
set(gca,'XScale','log')
xlim([.4,5.6])
ylim([0.5,15.5])
fill([.4 5.6 5.6 .4],[14.5 14.5 15.5 15.5],[.7 .7 .7],'FaceAlpha',0.2,'EdgeColor',[.7 .7 .7],'EdgeAlpha',0.2)
title('Systolic blood pressure')
xlabel('OR per 100mg daily use of sildenafil')
print(f,'./Figures/T2MR_SBP.png','-r300','-dpng')

%% ========================================================================
% Instruments effect size - harmonised data
% =========================================================================
clear all; close all;
names = {'DBP','SBP'};
names1 = {'Diastolic','Systolic'};
for i = 1:2
    t1 = readtable([names{i} '_AlzD_harmonised.txt']);
    t2 = readtable([names{i} '_AlzD_harmonised_2.txt']);
    
    blue = [114 172 185]./255;
    or   = [239 186 0]./255;

    f = figure(1); clf; hold on; grid on; box on;

    errorbar(flip(t1.Var11), 1:6-i, flip(-1.96*t1.Var13), flip(1.96*t1.Var13),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",'k',"Color",'k')
    errorbar(flip(t1.Var12), 1.05:6.05-i, flip(-1.96*t1.Var14), flip(1.96*t1.Var14),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
    errorbar(flip(t2.Var12), 0.95:5.95-i, flip(-1.96*t2.Var14), flip(1.96*t2.Var14),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or,"Color",or)
    xline(0,'LineStyle','--','LineWidth',2,'Color','k');
    ax = gca; ax.YTick = [1:6-i]; ax.YTickLabel = flip(t1.Var2);
    ylim([0.5,6.5-i])
    xlim([-.5,0.75])
    title([names1{i} ' blood pressure']);
    xlabel('Effect size [mmHg] / Effect size [log(OR)]')
    l = legend('IV-Blood pressure','IV-AlzD (Lambert)', 'IV-AlzD (Wightman)','Location','Northoutside','NumColumns',2);
    print(f,['./Figures/Instruments_' names{i} '.png'],'-r300','-dpng')
end


