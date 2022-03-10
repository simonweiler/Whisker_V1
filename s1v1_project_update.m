%% Load data structure for S1 V1
str   = 'D:\Postdoc_Margrie\Projects\Whisker\output_structure';
folder_list = uipickfiles('FilterSpec',str);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% Readout of cells 
%get all excitatory cells with EPSC/IPSC L2/3
temp1=[];
lv=[0 1];
for i=1:2
temp1(i,:)=cell_selecter(Ephys,'label',lv(i),'sol',2,'layer',3);
end
pyr_cs=sum(temp1);
%  get all cre on/ cre off L2/3 NON PAIRED
pyr_cs_creon=cell_selecter(Ephys,'label',1,'sol',2,'layer',3);
pyr_cs_creoff=cell_selecter(Ephys,'label',0,'sol',2,'layer',3);
% Cs-gluc cells in L23 PAIRED
temp1=[];temp2=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',2,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',2,'geno',7,'pair',i);
end
cre_on_cs=sum(temp1);cre_off_cs=sum(temp2);
%Paired Cs-gluc cells in L23
temp1=[];temp2=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',1,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',1,'geno',7,'pair',i);
end
cre_on_k=sum(temp1);cre_off_k=sum(temp2);
%% Plot example 
cnr=33 %210914SW0002 (nonlabelled)
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

%% Cre on vs cre off 
%CRE ON
cnr=7;%210907SW001
ov_min=-150;ov_max=20;
range=2000:9000;
temp=[];temp=find(pyr_cs_creon==1);
fig4=figure;set(fig4, 'Position', [200, 200, 250, 200]);set(gcf,'color','w');
subplot(1,2,1)
if max(Ephys(temp(cnr)).sub_traces_high(range,2))>max(Ephys(temp(cnr)).sub_traces_high(range,1))
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','m','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;
else
 plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','m','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
ylim([ov_min-10 ov_max]);title('Cre+','Color','m');

%CRE OFF
cnr=4;%210907SW001
temp=[];temp=find(pyr_cs_creoff==1);
subplot(1,2,2);
if max(Ephys(temp(cnr)).sub_traces_high(range,2))>max(Ephys(temp(cnr)).sub_traces_high(range,1))
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
 plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
ylim([ov_min-10 ov_max]);title('Cre-','Color','k');


%% Paired comparison EPSC IPSC using the first high frequcny pulse
[epsc_on_train ipsc_on_train e_i_ratio_on_train] = readout_amp(Ephys,cre_on_cs ,2);
[epsc_off_train ipsc_off_train e_i_ratio_off_train] = readout_amp(Ephys,cre_off_cs ,2);
[epsc_on_traink tr trtt] = readout_amp(Ephys,cre_on_k ,2);
[epsc_off_traink trr trtttt] = readout_amp(Ephys,cre_off_k ,2);
% only using pairs
cl={'m','k'};
data=[];data=[[epsc_on_train epsc_on_traink]' [epsc_off_train epsc_off_traink]'];
paired_plot_box(data,cl);
data=[];data=[[ipsc_on_train]' [ipsc_off_train]'];
cl={'m','k'};
paired_plot_box(data,cl);

%% Time to peak ex and in for retro cells
%
temp=[];
temp=find(pyr_cs==1);
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

%% Show E/I ratio from all PYR cells
[epsc_pyr_long ipsc_pyr_long e_i_ratio_pyr_long] = readout_amp(Ephys,pyr_cs ,1);
[epsc_pyr_hf ipsc_pyr_hf e_i_ratio_pyr_hf] = readout_amp(Ephys,pyr_cs ,2);
[epsc_pyr_hf2 ipsc_pyr_hf2 e_i_ratio_pyr_hf2] = readout_amp(Ephys,pyr_cs ,3);
 edges = [0:0.25:2];
 fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
 h1=histogram(e_i_ratio_pyr_long,edges,'Normalization','probability');h1.FaceColor='k'
hold on;h2=histogram(e_i_ratio_pyr_hf,edges,'Normalization','probability');h2.FaceColor='g'
 hold on;h3=histogram(e_i_ratio_pyr_hf2,edges,'Normalization','probability');h3.FaceColor='m'
 box off;xlabel('E/I ratio');ylabel('Relative counts');legend({'1 Hz','25 Hz','50 Hz'});
 hold on;plot([1 1],[0 0.3],'--k');legend boxoff  
 %% alternative; just show one with median 
  fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
hold on;h2=histogram(e_i_ratio_pyr_hf,edges,'Normalization','probability');h2.FaceColor='k';h2.FaceAlpha=0.5;
 box off;xlabel('E/I ratio');ylabel('Relative counts');
 hold on;plot([1 1],[0 0.4],'--k');
 hold on;plot([nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf))) nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf)))],[0.4 0.4],...
     'Marker','v','MarkerFaceColor','k','MarkerEdgeColor','k');title('25 Hz');set(gca,'FontSize',10);
%% Compare in bar plot all different type of stimuli 
a1=[];a1=find(~isinf(e_i_ratio_pyr_long)==1 & ~isnan(e_i_ratio_pyr_long)==1) ; 
a2=[];a2=find(~isinf(e_i_ratio_pyr_hf)==1 & ~isnan(e_i_ratio_pyr_hf)==1) ; 
a3=[];a3=find(~isinf(e_i_ratio_pyr_hf2)==1 & ~isnan(e_i_ratio_pyr_hf2)==1) ; 


fig6= figure;set(fig6, 'Name', 'compare EI');set(fig6, 'Position', [200, 300, 150, 250]);set(gcf,'color','w');
gr_m=[nanmean(e_i_ratio_pyr_long(a1)) nanmean(e_i_ratio_pyr_hf(a2))...
    nanmean(e_i_ratio_pyr_hf2(a3))];  
gr_sem=[nanstd(e_i_ratio_pyr_long(a1))/sqrt(length(e_i_ratio_pyr_long(a1)))...
    nanstd(e_i_ratio_pyr_hf(a2))/sqrt(length(e_i_ratio_pyr_hf(a2)))...
     nanstd(e_i_ratio_pyr_hf2(a3))/sqrt(length(e_i_ratio_pyr_hf2(a3)))];
hold on;
b=bar(1,gr_m(1));b.FaceColor='k';b.FaceAlpha=0.75;
b=bar(2,gr_m(2));b.FaceColor='k';b.FaceAlpha=0.5;
b=bar(3,gr_m(3));b.FaceColor='k';b.FaceAlpha=0.25;
% hold on;scatter(ones(length(e_i_ratio_pyr_long(a1)),1),e_i_ratio_pyr_long(a1),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf(a2)),1)*2,e_i_ratio_pyr_hf(a2),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf2(a3)),1)*3,e_i_ratio_pyr_hf2(a3),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
hold on;er=errorbar(1:3,gr_m,gr_sem);er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
text(0.7,0.1,num2str(length(e_i_ratio_pyr_long(a1))),'Color','w');hold on;
text(1.7,0.1,num2str(length(e_i_ratio_pyr_hf(a2))),'Color','w');hold on;
text(2.7,0.1,num2str(length(e_i_ratio_pyr_hf2(a3))),'Color','w');hold on;
xticks([1:1:3]);ylabel('E / I ratio');xticklabels({'1Hz','25Hz','50Hz'});xtickangle(45);set(gca,'FontSize',10);
%% compare input to PV and PYR
%CS solution 
pyr_cs_pv=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5);
pv_cs_pv=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5);
%K solution only pairs
temp1=[];temp2=[];
for i=1:3
temp1(i,:) = cell_selecter(Ephys,'label',[0],'sol',1,'geno',5,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[3],'sol',1,'geno',5,'pair',i);
end
pyr_k_pv=sum(temp1);
pv_k_pv=sum(temp2);

%Readout peaks 
[epsc_pyr ipsc_pyr e_i_ratio_pyr] = readout_amp(Ephys,pyr_cs_pv ,2);
[epsc_pv ipsc_pv e_i_ratio_pv] = readout_amp(Ephys,pv_cs_pv ,2);
[epsc_pyr_k] = readout_amp(Ephys,pyr_k_pv ,2);
[epsc_pv_k] = readout_amp(Ephys,pv_k_pv ,2);

%EPSC
data=[];data=[[epsc_pyr epsc_pyr_k]' [epsc_pv epsc_pv_k]'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PYR','PV'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);
%IPSC
data=[];data=[[ipsc_pyr]' [ipsc_pv]'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PYR','PV'});ylabel('IPSC Amplitude (pA)');set(gca,'FontSize',10);
%E/I ratio
data=[];data=[[e_i_ratio_pyr]' [e_i_ratio_pv]'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PYR','PV'});ylabel('E / I ratio');set(gca,'FontSize',10);
%% extract normalized traces EPSC PYR
range=[];range=1:60000;
temp=[];temp=find(pyr_cs_pv==1 | pyr_k_pv==1);
range_save=[];range_save=5000:6000;
traces_pyr=[];
for i=1:length(temp)
 if   size(Ephys(temp(i)).sub_traces_high(range,:),2)>1 
      
if isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,1)<-50))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,2))>max(Ephys(temp(i)).sub_traces_high(range,1))
traces_pyr(:,i)=Ephys(temp(i)).sub_traces_high(range_save,1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,1)));
elseif isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,2)<-50))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,1))>max(Ephys(temp(i)).sub_traces_high(range,2))
traces_pyr(:,i)=Ephys(temp(i)).sub_traces_high(range_save,2)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,2)));
else
    traces_pyr(:,i)=ones(length(range_save),1)*NaN;
end
 else
     if sum(Ephys(temp(i)).high_n)==0
         traces_pyr(:,i)=ones(length(range_save),1)*NaN;
     else
        traces_pyr(:,i)=Ephys(temp(i)).sub_traces_high(range_save,1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,1)));
     end
 end
end
%% extract normalized traces EPSC PV
range=[];range=1:60000;
temp=[];temp=find(pv_cs_pv==1 | pv_k_pv==1);
range_save=[];range_save=5000:6000;
traces_pv=[];
for i=1:length(temp)
 if   size(Ephys(temp(i)).sub_traces_high(range,:),2)>1 
      
if isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,1)<-50))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,2))>max(Ephys(temp(i)).sub_traces_high(range,1))
traces_pv(:,i)=Ephys(temp(i)).sub_traces_high(range_save,1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,1)));
elseif isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,2)<-50))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,1))>max(Ephys(temp(i)).sub_traces_high(range,2))
traces_pv(:,i)=Ephys(temp(i)).sub_traces_high(range_save,2)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,2)));
else
    traces_pv(:,i)=ones(length(range_save),1)*NaN;
end
 else
     if isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,1)<-50))==1
        traces_pv(:,i)=ones(length(range_save),1)*NaN;    
     else
         traces_pv(:,i)=Ephys(temp(i)).sub_traces_high(range_save,1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,1)));
     end
 end
end

%% show example trace 
%PYR
cnr=2;%
ov_min=-500;ov_max=50;
range=2000:9000;
temp=[];temp=find(pyr_k_pv==1);
fig4=figure;set(fig4, 'Position', [200, 200, 250, 200]);set(gcf,'color','w');
subplot(1,2,1)
%if max(Ephys(temp(cnr)).sub_traces_high(range,2))>max(Ephys(temp(cnr)).sub_traces_high(range,1))
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;
% else
%  plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','m','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('PYR','Color','k');

%PV
cnr=2;%210907SW001
temp=[];temp=find(pv_k_pv==1);
subplot(1,2,2);
%if max(Ephys(temp(cnr)).sub_traces_high(range,2))>max(Ephys(temp(cnr)).sub_traces_high(range,1))
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
% else
%  plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','k','LineWidth',1);set(gca,'box','off');
%hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
% end
ylim([ov_min-10 ov_max]);title('PV','Color',[0.8500 0.3250 0.0980]);
% 
%% Spiking or not PYR vs nonFS and FS
temp1=[];lab=[];lab=[0 1];
for i=1:2
temp1(i,:) = cell_selecter(Ephys,'label',lab(i),'sol',1);
end
pyr_k=sum(temp1);

temp2=[];lab=[];lab=[2 3];
for i=1:2
temp2(i,:) = cell_selecter(Ephys,'label',lab(i),'sol',1);
end
in_k=sum(temp2);
epsp_pyr=[];epsp_in=[];
[epsp_pyr] = readout_amp_epsp(Ephys,pyr_k ,2,sr);
[epsp_in] = readout_amp_epsp(Ephys,in_k ,2,sr);

temp=[];maxsF=[];
temp=find(in_k==1);
pv_label=[Ephys(temp).label]==3;
for i=1:length(temp)
    try
    maxsF(i)=max(Ephys(temp(i)).IV.spikecount);
    catch
    maxsF(i)=NaN;
    end
end
temp=[];maxsF_pyr=[];
temp=find(pyr_k==1);
for i=1:length(temp)
    try
    maxsF_pyr(i)=max(Ephys(temp(i)).IV.spikecount);
    catch
    maxsF_pyr(i)=NaN;
    end
end
%bug in first cell
maxsF(1)=86;
%fraction of spikers
spike_FS=sum(epsp_in>50)/(sum(maxsF>50)+1);
spike_nFS=0;
spike_PYR=nansum(epsp_pyr>50)/length((epsp_pyr));
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 150, 250]);set(gcf,'color','w');
gr_m=[spike_PYR spike_nFS spike_FS ];  
hold on;
b=bar(1,gr_m(1));b.FaceColor='k';b.FaceAlpha=0.75;
b=bar(2,gr_m(2));b.FaceColor=[0.8500 0.3250 0.0980];b.FaceAlpha=0.25;
b=bar(3,gr_m(3));b.FaceColor=[0.8500 0.3250 0.0980];b.FaceAlpha=0.75;
% text(0.7,0.1,num2str(length(e_i_ratio_pyr_long(a1))),'Color','w');hold on;
% text(1.7,0.1,num2str(length(e_i_ratio_pyr_hf(a2))),'Color','w');hold on;
% text(2.7,0.1,num2str(length(e_i_ratio_pyr_hf2(a3))),'Color','w');hold on;
xticks([1:1:3]);ylabel('Fraction of cells firing APs');xticklabels({'PYR','nFS IN','FS IN'});xtickangle(45);set(gca,'FontSize',10);
ylim([0 1]);

% scatter plot
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 230, 230]);set(gcf,'color','w');
scatter(maxsF_pyr,epsp_pyr,25,'o','MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
scatter(maxsF,epsp_in,25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980]);
ylabel('EPSP amplitude (mV)');xlabel('max Spike frequency (Hz)');legend({'PYR','IN'}); set(gca,'FontSize',10);legend boxoff
%% showing EPSPs
a=[];a=find(maxsF>50 | pv_label==1);
b=[];b=find(maxsF<50);
data=[];data=[epsp_pyr epsp_in(b) epsp_in(a)]';
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
 hold on;scatter(ones(length(epsp_pyr),1),epsp_pyr,25,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
 hold on;scatter(ones(length(epsp_in(b)),1)*2,epsp_in(b),25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5]);
 hold on;scatter(ones(length(epsp_in(a)),1)*3,epsp_in(a),25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980]);
 hold on;plot([0.6 1.4],[nanmean(epsp_pyr) nanmean(epsp_pyr)],'Color','k');
 hold on;plot([1.6 2.4],[nanmean(epsp_in(b)) nanmean(epsp_in(b))],'Color',[0.5 0.5 0.5]);
 hold on;plot([2.6 3.4],[nanmean(epsp_in(a)) nanmean(epsp_in(a))],'Color',[0.8500 0.3250 0.0980]);
  %hold on;p=plot([0 4],[50 50],':k');
 xlim([0 4]);xticks([1:1:3]);ylabel('EPSP amplitude (mV)');xticklabels({'PYR','nFS IN','FS IN'});xtickangle(45);
 set(gca,'FontSize',10);
%% exmaple EPSP traces
cnr=9;%
ov_min=-5;ov_max=100;
range=2000:9000;
temp=[];temp=find(pyr_k==1);
fig4=figure;set(fig4, 'Position', [200, 200, 150, 300]);set(gcf,'color','w');
subplot(3,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','k','LineWidth',1);set(gca,'box','off');
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('PYR','Color','k');

cnr=3;%
ov_min=-5;ov_max=100;
range=2000:9000;
temp=[];temp=find(in_k==1);
subplot(3,1,2)
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color',[0.5 0.5 0.5],'LineWidth',1);set(gca,'box','off');
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('nFS IN','Color',[0.5 0.5 0.5]);

cnr=8;%
ov_min=-5;ov_max=100;
range=2000:9000;
temp=[];temp=find(in_k==1);
subplot(3,1,3)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1);set(gca,'box','off');
hold on;plot([0.15*sr 0.15*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('FS IN','Color',[0.8500 0.3250 0.0980]);

%% Intrinsic properties
[rmp_pyr maxsp_pyr rheo_pyr rin_pyr tau_pyr sag_pyr trace_pyr spike_time_pyr] = passive_readout(Ephys,pyr_k);
[rmp_in maxsp_in rheo_in rin_in tau_in sag_in trace_in spike_time_in] = passive_readout(Ephys,in_k);
%% 
p1=rin_pyr;p2=rin_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
ylabel('Input resistance');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

p1=rheo_pyr;p2=rheo_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
ylabel('Threshold current (pA)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

p1=rmp_pyr;p2=rmp_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
ylabel('RMP (mV)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

%% Active parameters
active_pyr=sp_parameters_pandora(trace_pyr,2);
active_in=sp_parameters_pandora(trace_in,2);
%% 
p1=[];p2=[];p1=active_pyr(10,:);p2=active_in(10,:);
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
ylabel('APVslope (mV/ms)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);
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
% %% 
