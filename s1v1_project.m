%% Load data structure for S1 V1
str   = 'D:\Postdoc_Margrie\Projects\Whisker\output_structure';
folder_list = uipickfiles('FilterSpec',str);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% get all excitatory cells with EPSC/IPSC L2/3
temp1=[];
lv=[0 1];
for i=1:2
temp1(i,:)=cell_selecter(Ephys,'label',lv(i),'sol',2,'layer',3);
end
pyr_cs=sum(temp1);
%% Plot example 
cnr=31 %210914SW0002 (nonlabelled)
ov_min=-400;ov_max=600;
temp=[];temp=find(pyr_cs==1);
fig4=figure;set(fig4, 'Position', [200, 200, 200, 300]);set(gcf,'color','w');
if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
%%  get all cre on/ cre off L2/3 NON PAIRED
pyr_cs_creon=cell_selecter(Ephys,'label',1,'sol',2,'layer',3);
pyr_cs_creoff=cell_selecter(Ephys,'label',0,'sol',2,'layer',3);
%% Cre on vs cre off 
%CRE ON
cnr=7;%210907SW001
ov_min=-150;ov_max=20;
temp=[];temp=find(pyr_cs_creon==1);
fig4=figure;set(fig4, 'Position', [200, 200, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;
else
 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
ylim([ov_min-10 ov_max]);title('Cre+','Color','r');

%CRE OFF
cnr=4;%210907SW001
temp=[];temp=find(pyr_cs_creoff==1);
subplot(1,2,2);
if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
ylim([ov_min-10 ov_max]);title('Cre-','Color','k');

%% Cs-gluc cells in L23 PAIRED
temp1=[];temp2=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',2,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',2,'geno',7,'pair',i);
end
cre_on_cs=sum(temp1);
cre_off_cs=sum(temp2);
%% 
%Paired Cs-gluc cells in L23
temp1=[];temp2=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',1,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',1,'geno',7,'pair',i);
end
cre_on_k=sum(temp1);
cre_off_k=sum(temp2);
%% 


%% use high frequency pulse for this analysis 
[epsc_on_train ipsc_on_train e_i_ratio_on_train] = readout_amp(Ephys,cre_on_cs ,2);
[epsc_off_train ipsc_off_train e_i_ratio_off_train] = readout_amp(Ephys,cre_off_cs ,2);
[epsc_on_traink tr trtt] = readout_amp(Ephys,cre_on_k ,2);
[epsc_off_traink trr trtttt] = readout_amp(Ephys,cre_off_k ,2);

% %% 
% edges = [0:0.25:2];
% fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
% h1=histogram(all_eit,edges,'Normalization','probability');h1.FaceColor='k'
% hold on;h2=histogram(all_eif,edges,'Normalization','probability');h2.FaceColor='g'
% hold on;h3=histogram(all_eif2,edges,'Normalization','probability');h3.FaceColor='m'
% box off;xlabel('E/I ratio');ylabel('Relative counts');legend({'1 Hz','25 Hz','50 Hz'});
% hold on;plot([1 1],[0 0.3],'--k');legend boxoff 
%% only using pairs
data=[];data=[[epsc_on_train epsc_on_traink]' [epsc_off_train epsc_off_traink]'];
paired_plot_box(data);
data=[];data=[[ipsc_on_train]' [ipsc_off_train]'];
paired_plot_box(data);
%% traces for epsc and ipsc for CRE ON and CRE OFF
%epsc
temp=[];
range=6;
temp=find(cre_on_cs);
for cnr=1:length(temp)
if max(Ephys(temp(cnr)).sub_traces_high(1:range*sr,2))>max(Ephys(temp(cnr)).sub_traces_high(1:range*sr,1))
     traces_epsc_creon(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,1);
     traces_ipsc_creon(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,2);
else
     traces_epsc_creon(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,2);
     traces_ipsc_creon(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,1);
end
end
%ipsc 
temp=[];
range=6;
temp=find(cre_off_cs);
for cnr=1:length(temp)
if max(Ephys(temp(cnr)).sub_traces_high(1:range*sr,2))>max(Ephys(temp(cnr)).sub_traces_high(1:range*sr,1))
     traces_epsc_creoff(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,1);
     traces_ipsc_creoff(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,2);
else
     traces_epsc_creoff(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,2);
     traces_ipsc_creoff(:,cnr)=Ephys(temp(cnr)).sub_traces_high(1:range*sr,1);
end
end

%% 



%% Plotexmaple epsc/ipsc single pulse
%cnr=22;
cnr=22
ov_min=-400;ov_max=600;
temp=[];
temp=find(all_cs==1);
fig4=figure;set(fig4, 'Position', [200, 200, 700, 300]);set(gcf,'color','w');
subplot(1,3,1)
if max(Ephys(temp(cnr)).sub_traces_train(1:5*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:5*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:5*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
subplot(1,3,2)
temp=[];
temp=find(all_cs==1);
if max(Ephys(temp(cnr)).sub_traces_high(1:3*sr,2))>max(Ephys(temp(cnr)).sub_traces_high(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_high(1:3*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_high(1:3*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_high(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_high(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end

subplot(1,3,3)
temp=[];
temp=find(all_cs==1);
if max(Ephys(temp(cnr)).sub_traces_highf(1:3*sr,2))>max(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_highf(1:3*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:3*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end

%% Single first pulse
cnr=22
ov_min=-400;ov_max=600;
temp=[];
temp=find(all_cs==1);
fig4=figure;set(fig4, 'Position', [200, 200, 200, 300]);set(gcf,'color','w');
if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
%% 
%% Time to peak ex and in for retro cells
%Ntsr1 mouse line, Cs-gluc, retro cells

temp=[];
temp=find(all_cs==1);
close all
for i=1:length(temp)
    if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
[t_ex(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:30:5000],[5000:30:6000],20,0);
[t_in(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:30:5000],[5000:30:6000],20,1);
    else
      [t_ex(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:30:5000],[5000:30:6000],20,0);
      [t_in(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:30:5000],[5000:30:6000],20,1);
    end
end
%% Plot EX vs IN delay difference
t_ein=[];
t_ein=[t_ex' t_in'];
[gh gm]=find(isnan(t_ein));
t_ein(unique(gh),:)=[];

cl={'r','b'};
data=[];data=t_ein;
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
%hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','r','Markersize',7);
hold on;errorbar([0.75],nanmean(data(:,1)),nanstd(data(:,1),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor','r','Markersize',7);
hold on;errorbar([2.25],nanmean(data(:,2)),nanstd(data(:,2),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor','b','Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'EX','IN'});ylabel('Delay (ms)');set(gca,'FontSize',10);
%% TTX wash in 
all_cs_ttx = cell_selecter(Ephys,'sol',2,'drugs',1);
temp=[];
temp=find(all_cs_ttx==1);
cnr=3
ov_min=-20;ov_max=300;
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,3),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,2),'Color','b','LineWidth',1.2);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
subplot(1,2,2)
plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,4),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,1),'Color','r','LineWidth',1.2);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%% IPSC quantification
cl={'b',[0.5 0.5 0.5]};
ttx_ipsc=[max(abs(Ephys(temp(3)).highf_p(:,2))) max(abs(Ephys(temp(3)).highf_p(:,3)));...
    max(abs(Ephys(temp(2)).highf_p(:,2))) max(abs(Ephys(temp(2)).highf_p(:,3)));...
  max(abs(Ephys(temp(1)).highf_p(:,1))) max(abs(Ephys(temp(1)).highf_p(:,3))) ];
data=[];data=ttx_ipsc;
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
end
hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
box off;set(gca,'FontSize',10);
 [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'no TTX','TTX'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
%% EPSC
ttx_epsc=[max(abs(Ephys(temp(3)).highf_n(:,1))) max(abs(Ephys(temp(3)).highf_n(:,4)));...
    max(abs(Ephys(temp(2)).highf_n(:,1))) max(abs(Ephys(temp(2)).highf_n(:,4)));...
  max(abs(Ephys(temp(1)).highf_n(:,3))) max(abs(Ephys(temp(1)).highf_n(:,4))) ];
data=[];data=ttx_epsc;
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
end
hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
box off;set(gca,'FontSize',10);
 [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'no TTX','TTX'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
%% TTX modulation index
par=[(ttx_ipsc(:,2)-ttx_ipsc(:,1))./(ttx_ipsc(:,2)+ttx_ipsc(:,1)); (ttx_epsc(:,2)-ttx_epsc(:,1))./(ttx_epsc(:,2)+ttx_epsc(:,1))]
s1=[1:3];s2=[4:6]
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'IPSC' ,'EPSC'});ylabel('TTX modulation index');set(gca,'FontSize',10);xtickangle(45);
%% 
all_red=[];all_nred=[];
all_red = cell_selecter(Ephys,'label',1,'pair',1,'geno',7);
all_nred = cell_selecter(Ephys,'label',0,'pair',1,'geno',7);
%% 
all_red=[];all_nred=[];
all_red = cell_selecter(Ephys,'label',1,'geno',7);
all_nred = cell_selecter(Ephys,'label',0,'geno',7);
%% 
red_epsc=[];red_epscf=[];nred_epsc=[];nred_epscf=[];
temp=[];
for i=1:length(find(all_red==1));
   temp=find(all_red==1);
   red_epsc(i)=max(abs(Ephys(temp(i)).train_n(:)));
end
temp=[];
   for i=1:length(find(all_nred==1));
   temp=find(all_nred==1);
   nred_epsc(i)=max(abs(Ephys(temp(i)).train_n(:)));
   end
   
temp=[];
for i=1:length(find(all_red==1));
   temp=find(all_red==1);
   red_epscf(i)=max(abs(Ephys(temp(i)).high_n(:)));
end
temp=[];
   for i=1:length(find(all_nred==1));
   temp=find(all_nred==1);
   nred_epscf(i)=max(abs(Ephys(temp(i)).high_n(:)));
   end
   
%    temp=[];
% for i=1:length(find(all_red==1));
%    temp=find(all_red==1);
%    red_epschf(i)=max(abs(Ephys(temp(i)).highf_n(:)));
% end
% temp=[];
%    for i=1:length(find(all_nred==1));
%    temp=find(all_nred==1);
%    nred_epschf(i)=max(abs(Ephys(temp(i)).highf_n(:)));
%    end
   
   %% All against all
   par=[red_epsc nred_epsc]
   s1=1:length(red_epsc)
   s2=length(red_epsc)+1:length(nred_epsc)+length(red_epsc)
  [statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'Cre+' ,'Cre-'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);xtickangle(45);

%% only using pairs
cl={'r',[0.5 0.5 0.5]};
data=[];data=[red_epscf' nred_epscf'];
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
%hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
%hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','r','Markersize',7);
hold on;errorbar([0.75],nanmean(data(:,1)),nanstd(data(:,1),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor','r','Markersize',7);
hold on;errorbar([2.25],nanmean(data(:,2)),nanstd(data(:,2),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor',[0.5 0.5 0.5],'Markersize',7);
 [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
xticklabels({'Cre+','Cre-'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);

%% 

cnr=6
ov_min=-400;ov_max=600;
temp=[];
temp=find(all_red==1);
fig4=figure;set(fig4, 'Position', [200, 200, 600, 300]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_high(1:2*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
ylim([-200 20])
ov_min=-400;ov_max=600;
title('Cre+')
temp=[];
temp=find(all_nred==1);
subplot(1,2,2)

plot(Ephys(temp(cnr)).sub_traces_high(1:2*sr,1),'Color',[0.5 0.5 0.5],'LineWidth',1);set(gca,'box','off');
ylim([-200 20])
title('Cre-')
%% 

 plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
