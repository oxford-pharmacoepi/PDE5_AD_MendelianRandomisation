clear all; close all;
t = readtable('mr_results.csv');

f = figure(1); hold on; grid on; box on;
errorbar(t.OR,1:2, t.OR-t.lower, t.upper-t.OR,"s","horizontal","LineWidth",1.25,"MarkerFaceColor",'k',"Color",'k')
fill([-1.5 2.5 1.2 -0.6],[12 12 10.5 10.5],[0.7 0.7 0.7],"FaceAlpha",0.2,"EdgeColor",[0.7 0.7 0.7],"EdgeAlpha",0.2)