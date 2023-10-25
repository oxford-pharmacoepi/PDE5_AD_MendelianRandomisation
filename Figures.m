% =========================================================================
% All the results
% =========================================================================
clear all; close all;
% First analysis
t1 = readtable('scaled_mr_results.csv');
% de Rojas
t2 = readtable('scaled_mr_results_3.csv');
% Wigthman
t3 = readtable('scaled_mr_results_2.csv');

blue = [114 172 185]./255;
or  = [213 113 0]./255;
or2  = [239 186 0]./255;

f = figure(1); hold on; grid on; box on;
errorbar(flip(t1.OR), [1.05, 1.45], flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), [1, 1.4], flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or,"Color",or)
errorbar(flip(t3.OR), [0.95, 1.35], flip(t3.OR-t3.lower), flip(t3.upper-t3.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or2,"Color",or2)

xline(1,'LineStyle','--','LineWidth',2)
ylim([.8,1.6])
ax = gca; ax.YTick = [1, 1.4]; ax.YTickLabel = {'SBP','DBP'}; ax.XLim = [0.4 5.5]; ax.XTick = [0.5,1,2,3,4,5]
l = legend('Lambert et al. 2013', 'de Rojas et al. 2021','Wightman et al. 2021','Location','Northoutside','NumColumns',3);
l.Position = [0.1    0.95    0.85   0.0464];
set(gca, 'XScale', 'log')
xlabel('OR per 100mg daily use of sildenafil')
print(f,'./Figures/All.png','-r300','-dpng')

%% ========================================================================
% Two step results // Diastolic blood pressure
% =========================================================================
clear all; close all;
a = 3

if a == 1
    % First analysis - main
    t1 = readtable('scaled_mr_results.csv');
    t1 = t1(1,:);
    t2 = readtable('scaled_mr_results_twoStep_DBP.csv');
    t2 = t2(1:14,:);
    blue = [114 172 185]./255;
    name = 'Lambert';
elseif a == 2
    t1 = readtable('scaled_mr_results_3.csv');
    t1 = t1(1,:);
    t2 = readtable('scaled_mr_results_twoStep_DBP_3.xlsx');
    blue  = [213 113 0]./255;
    name = 'deRojas';
elseif a == 3
    t1 = readtable('scaled_mr_results_2.csv');
    t1 = t1(1,:);
    t2 = readtable('scaled_mr_results_twoStep_DBP_2.csv');
    blue  = [239 186 0]./255;
    name = 'Wightman';
end

f = figure(1); hold on; grid on; box on;
xline(1,'LineStyle','--','LineWidth',2,'Color','k')
errorbar(flip(t1.OR), 15, flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), 1:14, flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
ax = gca; ax.YTick = [1:15]; ax.YTickLabel = flip({'Main','BMI','Impedance of leg (right)','Impedance of leg (left)','Impedance of arm (right)','Impedance of arm (left)','Impedance of whole body','Standing height','Plateletcrit','Myeloid white cell count','Platelet count','White blood cell count','Coronary artery disease','Granulocyte count','Sum basophil neutrophil counts'});
set(gca,'XScale','log')
if a == 1
    ax.XTick = [.5 .75 1 1.5 2];
    xlim([.45,2.3])
elseif a == 2
    ax.XTick = [.5 .75 1 1.5 2];
    xlim([.45,2.3])
elseif a == 3
    ax.XTick = [.75 1 1.5 2 3 4];
    xlim([.7,4.25])
end
 fill([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[14.5 14.5 15.5 15.5],[.7 .7 .7],'FaceAlpha',0.2,'EdgeColor',[.7 .7 .7],'EdgeAlpha',0.2)
ylim([0.5,15.5])
title('Diastolic blood pressure')
xlabel('OR per 100mg daily use of sildenafil')
print(f,['./Figures/T2MR_DBP' name '.png'],'-r300','-dpng')

%% ========================================================================
% Systolic blood pressure
% =========================================================================
clear all; close all;
a = 3
if a == 1
    % First analysis - main
    t1 = readtable('scaled_mr_results.csv');
    t1 = t1(2,:);
    t2 = readtable('scaled_mr_results_twoStep_SBP.csv');
    t2 = t2(15:end,:);
    blue = [114 172 185]./255;
    name = 'Lambert';
elseif a == 2
    t1 = readtable('scaled_mr_results_3.csv');
    t1 = t1(2,:);
    t2 = readtable('scaled_mr_results_twoStep_SBP_3.xlsx');
    blue  = [213 113 0]./255;
    name = 'deRojas';
elseif a == 3
    t1 = readtable('scaled_mr_results_2.csv');
    t1 = t1(2,:);
    t2 = readtable('scaled_mr_results_twoStep_SBP_2.csv');
    blue  = [239 186 0]./255;
    name = 'Wightman';
end

f = figure(1); hold on; grid on; box on;
xline(1,'LineStyle','--','LineWidth',2,'Color','k')
errorbar(flip(t1.OR), 15, flip(t1.OR-t1.lower), flip(t1.upper-t1.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
errorbar(flip(t2.OR), 1:14, flip(t2.OR-t2.lower), flip(t2.upper-t2.OR),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
ax = gca; ax.YTick = [1:15]; ax.YTickLabel = flip({'Main','BMI','Impedance of leg (right)','Impedance of leg (left)','Impedance of arm (right)','Impedance of arm (left)','Impedance of whole body','Standing height','Plateletcrit','Myeloid white cell count','Platelet count','White blood cell count','Coronary artery disease','Granulocyte count','Sum basophil neutrophil counts'});
set(gca,'XScale','log')
if a == 1
    ax.XTick = [.5 1 2 3 4 5];
    xlim([.4,5.6])
elseif a == 2
    ax.XTick = [.5 .75 1 1.5 2];
    xlim([.45,2.3])
elseif a == 3
    ax.XTick = [.75 1 1.5 2 3 4];
    xlim([.7,4.25])
end
ylim([0.5,15.5])
fill([ax.XLim(1) ax.XLim(2) ax.XLim(2) ax.XLim(1)],[14.5 14.5 15.5 15.5],[.7 .7 .7],'FaceAlpha',0.2,'EdgeColor',[.7 .7 .7],'EdgeAlpha',0.2)
title('Systolic blood pressure')
xlabel('OR per 100mg daily use of sildenafil')
print(f,['./Figures/T2MR_SBP_' name '.png'],'-r300','-dpng')

%% ========================================================================
% Instruments effect size - harmonised data
% =========================================================================
clear all; close all;
names = {'DBP','SBP'};
names1 = {'Diastolic','Systolic'};
for i = 1:2
    t1 = readtable([names{i} '_AlzD_harmonised.txt']);
    t3 = readtable([names{i} '_AlzD_harmonised_3.txt']);
    t2 = readtable([names{i} '_AlzD_harmonised_2.txt']);

    blue = [114 172 185]./255;
    or  = [213 113 0]./255;
    or2  = [239 186 0]./255;

    f = figure(1); clf; hold on; grid on; box on;

    errorbar(flip(t1.Var11), 1:6-i, flip(-1.96*t1.Var13), flip(1.96*t1.Var13),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",'k',"Color",'k')
    errorbar(flip(t1.Var12), 1.15:6.15-i, flip(-1.96*t1.Var14), flip(1.96*t1.Var14),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",blue,"Color",blue)
    errorbar(flip(t3.Var9), 1:6-i, flip(-1.96*t3.Var11), flip(1.96*t3.Var11),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or,"Color",or)
    errorbar(flip(t2.Var12), 0.85:5.85-i, flip(-1.96*t2.Var14), flip(1.96*t2.Var14),"s","horizontal","LineWidth",1.25,"MarkerFaceColor",or2,"Color",or2)

    xline(0,'LineStyle','--','LineWidth',2,'Color','k');
    ax = gca; ax.YTick = [1:6-i]; ax.YTickLabel = flip(t1.Var2);
    ylim([0.5,6.5-i])
    xlim([-.5,0.75])
    title([names1{i} ' blood pressure']);
    xlabel('Effect size [mmHg] / Effect size [log(OR)]')
    l = legend('IV-Blood pressure','IV-AlzD (Lambert)', 'IV-AlzD (de Rojas)','IV-AlzD (Wightman)','Location','Northoutside','NumColumns',2);
    print(f,['./Figures/Instruments_' names{i} '.png'],'-r300','-dpng')
end


