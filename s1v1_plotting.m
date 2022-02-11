%% Load data structure for S1V1
str_S1V1    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_S1V1);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;

%% use filter function 'cell selecter' to read out desired cells/line etc.
%Ntsr1 mouse line, K-gluc, retro cells
antero_cells = cell_selecter(Ephys, 'label',1, 'geno',7);
non_antero_cells = cell_selecter(Ephys, 'label',0, 'geno',7);
%% Plot example IV
cnr=7;step=[1:2:15];
temp=[];
temp=find(antero_cells==1);
ov_min=-150;ov_max=50;
fig4=figure;set(fig4, 'Position', [200, 800, 350, 300]);set(gcf,'color','w');
%plot(Ephys(temp(cnr)).IV.traces(:,step),'Color','k','LineWidth',1);set(gca,'box','off');
plot(Ephys(temp(cnr)).Ramp.traces,'Color','k','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-300*srF,Ephys(temp(cnr)).IV.RMP,[num2str(Ephys(temp(cnr)).IV.RMP),'mV'],'FontSize',9);
hold on;title('Antero+');axis off;
%% 
cnr=2;step=[1:2:15];
temp=[];
temp=find(non_antero_cells==1);
ov_min=-150;ov_max=50;
fig4=figure;set(fig4, 'Position', [200, 800, 350, 300]);set(gcf,'color','w');
%plot(Ephys(temp(cnr)).IV.traces(:,step),'Color','k','LineWidth',1);set(gca,'box','off');
plot(Ephys(temp(cnr)).Ramp.traces,'Color','k','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-300*srF,Ephys(temp(cnr)).IV.RMP,[num2str(Ephys(temp(cnr)).IV.RMP),'mV'],'FontSize',9);
hold on;title('Antero-');axis off;

%% Plot example input trace VC
cnr=8;
ov_min=-100;ov_max=2;
temp=[];
temp=find(antero_cells==1);
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');

plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','k','LineWidth',1);set(gca,'box','off');
%plot(Ephys(temp(cnr)).sub_traces_high(1:1*sr,1),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('Antero+');;ylabel('EPSC (mV)');xlabel('Time (s)');xticks([0:5000:20000]);xticklabels({'0','0.25','0.5','0.75','1'})
% %Scale bar
%  scale_x= 200;
%  scale_y= 40;
% %  %scale barx
%   hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
% %  %scale bary
%   hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
%% 

cnr=2;
hold on;
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
ov_min=-100;ov_max=2;
temp=[];
temp=find(non_antero_cells==1);
%plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,2),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');
hold on;;title('Antero-','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylabel('EPSP (mV)');xlabel('Time (s)');xticklabels({'0','0.25','0.5','0.75','1'})
% axis off;
% %Scale bar
%  scale_x= 200;
%  scale_y= 40;
%  %scale barx
%  hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
%  %scale bary
%  hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
%% 

cnr=10;
ov_min=-200;ov_max=300;
temp=[];
temp=find(antero_cells==1);
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');

% plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
 plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,end),'Color','r','LineWidth',1);set(gca,'box','off');
 hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,end-1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('Antero+');;ylabel('PSC (pA)');xlabel('Time (s)');xticks([0:5000:20000]);xticklabels({'0','0.25','0.5','0.75','1'})
% %Scale bar
%  scale_x= 200;
%  scale_y= 40;
% %  %scale barx
%   hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
% %  %scale bary
%   hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
%% 

cnr=3;
ov_min=-0.5;ov_max=2;
temp=[];
temp=find(non_antero_cells==1);
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');

plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,5),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,5),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min ov_max]);
%  plot(Ephys(temp(cnr)).sub_traces_high(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
%  hold on;plot(Ephys(temp(cnr)).sub_traces_high(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('Antero-');;ylabel('PSC (pA)');xlabel('Time (s)');xticks([0:5000:20000]);xticklabels({'0','0.25','0.5','0.75','1'})
% %Scale bar
%  scale_x= 200;
%  scale_y= 40;
% %  %scale barx
%   hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
% %  %scale bary
%   hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
