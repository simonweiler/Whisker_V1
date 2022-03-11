%% Load data structure for S1 V1 using uipickfiles 
str   = 'D:\Postdoc_Margrie\Projects\Whisker\output_structure';
folder_list = uipickfiles('FilterSpec',str);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% Readout indices of different cells using cell selecter function 

%get all excitatory cells with EPSC/IPSC L2/3 with no TTX 4AP
temp1=[];temp2=[];pyr_cs=[];
lv=[0 1];
for i=1:2
temp1(i,:)=cell_selecter(Ephys,'label',lv(i),'sol',2,'layer',3,'drugs',0);
temp2(i,:)=cell_selecter(Ephys,'label',lv(i),'sol',2,'layer',3,'drugs',1);
end
pyr_cs=sum([temp1; temp2]);


% get all cre on/ cre off L2/3 NON PAIRED no TTX 4AP also including K and
% Cs
temp1=[];temp2=[];pyr_cs_creoff=[];pyr_cs_creon=[];pyr_k_creon=[];pyr_k_creoff=[];
temp1=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'drugs',0);
temp2=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'drugs',1);
pyr_cs_creoff=sum([temp1; temp2]);
pyr_cs_creon=cell_selecter(Ephys,'label',1,'sol',2,'layer',3);
pyr_k_creon=cell_selecter(Ephys,'label',1,'layer',3,'geno',7,'sol',1);
pyr_k_creoff=cell_selecter(Ephys,'label',0,'layer',3,'geno',7,'sol',1);


% Cs-gluc cells in L23 PAIRED
temp1=[];temp2=[];cre_on_cs=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',2,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',2,'geno',7,'pair',i);
end
cre_on_cs=sum(temp1);cre_off_cs=sum(temp2);


%K-gluc cells in L23 PAIRED
temp1=[];temp2=[];cre_on_k=[];
for i=1:5
temp1(i,:) = cell_selecter(Ephys,'label',[1],'sol',1,'geno',7,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',[0],'sol',1,'geno',7,'pair',i);
end
cre_on_k=sum(temp1);cre_off_k=sum(temp2);
%% PN vs IN (GAD/PV) read out using cell selecter function 

%NON PAIRED CS solution without TTX 4AP, including wash in
temp1=[];temp2=[];pyr_cs_pv_all=[];
temp1=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',0);
temp2=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',1);
pyr_cs_pv_all=sum([temp1; temp2]);
temp1=[];temp2=[];pv_cs_pv_all=[];
temp1=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',0);
temp2=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',1);
pv_cs_pv_all=sum([temp1; temp2]);


%PAIRED CS solution without TTX 4AP, excluding wash in 
pyr_cs_pv=[];pv_cs_pv=[];
pyr_cs_pv=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',0);
pv_cs_pv=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',0);


%PAIRED CS solution with TTX 4AP, be careful to place if statment to read traces
%out after wash in 
temp1=[];temp2=[];pyr_ttx=[];
temp1=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',1);
temp2=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',2);
pyr_ttx=sum([temp1; temp2]);
temp1=[];temp2=[];pv_ttx=[];
temp1=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',1);
temp2=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',2);
pv_ttx=sum([temp1; temp2]);

%PAIRED CS solution without TTX 4AP, excluding wash in 
pyr_cs_pv_washin=[];pv_cs_pv_washin=[];
pyr_cs_pv_washin=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'geno',5,'drugs',1);
pv_cs_pv_washin=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'geno',5,'drugs',1);

%TTX 4AP WASH IN experiments 
temp1=[];temp2=[];pyr_washin=[];
temp1=cell_selecter(Ephys,'label',0,'sol',2,'layer',3,'drugs',1);
temp2=cell_selecter(Ephys,'label',1,'sol',2,'layer',3,'drugs',1);
pyr_washin=sum([temp1; temp2]);
temp1=[];temp2=[];pv_washin=[];
temp1=cell_selecter(Ephys,'label',3,'sol',2,'layer',3,'drugs',1);
pv_washin=sum([temp1; temp2]);


%PAIRED K solution without TTX 4AP
temp1=[];temp2=[];pyr_k_pv=[];pv_k_pv=[];
for i=1:3
temp1(i,:) = cell_selecter(Ephys,'label',0,'sol',1,'geno',5,'pair',i);
temp2(i,:) = cell_selecter(Ephys,'label',3,'sol',1,'geno',5,'pair',i);
end
pyr_k_pv=sum(temp1);
pv_k_pv=sum(temp2);




%Spiking PN (cre off + cre on) vs IN (GAD+PV lines)
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

%%  Read put peaks from train for PN cs for train (long), middle frequency (high), highest (highf)
[epsc_pyr_long ipsc_pyr_long e_i_ratio_pyr_long] = readout_amp(Ephys,pyr_cs ,1,2);
[epsc_pyr_hf ipsc_pyr_hf e_i_ratio_pyr_hf] = readout_amp(Ephys,pyr_cs ,2,2);
[epsc_pyr_hf2 ipsc_pyr_hf2 e_i_ratio_pyr_hf2] = readout_amp(Ephys,pyr_cs ,3,2);
%% Time to peak ex and in for retro cells using first pulse of long train 
temp=[];t_ex=[];t_in=[];t_in_ex=[];trace_smooth_ex=[];trace_smooth_in=[];
temp=find(pyr_cs==1);
fc_1=3.5;
close all
for i=1:length(temp)
    try
    if max(Ephys(temp(i)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(i)).sub_traces_train(1:1*sr,1))==1
[t_ex(i) trace_smooth_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,0,fc_1);
[t_in(i) trace_smooth_in(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,1,fc_1);
[t_in_ex(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,0,fc_1);
    else
      [t_ex(i) trace_smooth_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,0,fc_1);
      [t_in(i) trace_smooth_in(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,1,fc_1);
      [t_in_ex(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,0,fc_1);
    end
    catch
     t_ex(i)=NaN;
     t_in(i)=NaN; 
     t_in_ex(i)=NaN;   
    end
end
close all;
%% weird cells : in train EPSC traces shows longer delay then at 0 mV because of thresholding, just missing the EPSC in the 0 mV trace
% fc_2;
% [t_in_ex(30)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,0,fc_2);
t_in_ex(30)=NaN;
t_ex(30)=NaN;
t_in(34)=NaN;
%% Main Figure panel c plot EPSC and ISPC example  
cnr=23 %210914SW0002 (nonlabelled cs solution, no drugs)
ov_min=-400;ov_max=600;temp=[];temp=find(pyr_cs==1);
fig4=figure;set(fig4, 'Position', [200, 200, 200, 300]);set(gcf,'color','w');
if max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1))
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
else
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
end
%% Pairecd compariosn of onset latency between EPSC and IPSC  panel d
t_ex_sub=[];t_ex_sub=t_ex(find(e_i_pyr_train>0));
t_in_sub=[];t_in_sub=t_in(find(e_i_pyr_train>0));
t_in_ex_sub=[];t_in_ex_sub=t_in_ex(find(e_i_pyr_train>0));
subex=[];subex=find(t_in_ex_sub<t_ex_sub);
t_ex_sub(subex)=t_in_ex_sub(subex);
t_ein=[];
t_ein=[t_ex_sub' t_in_sub'];
[gh gm]=find(isnan(t_ein));
t_ein(unique(gh),:)=[];
cl={'r','b'};
data=[];data=t_ein;;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);

%% %% TTX modulation index for PN ; not happy with readout
temp=[];temp=find(pyr_washin==1);ttx_ipsc=[];ttx_epsc=[];

ttx_ipsc=[max(abs(Ephys(temp(1)).high_p(1:2,2))) max(abs(Ephys(temp(1)).high_p(1:2,3)));...
    max(abs(Ephys(temp(2)).high_p(1:2,2))) max(abs(Ephys(temp(2)).high_p(1:2,3)));...
  max(abs(Ephys(temp(3)).high_p(1:2,2))) max(abs(Ephys(temp(3)).high_p(1:2,3)));...
   max(abs(Ephys(temp(4)).high_p(1:2,2))) max(abs(Ephys(temp(4)).high_p(1:2,4)))];

ttx_epsc=[max(abs(Ephys(temp(1)).high_n(1:2,1))) max(abs(Ephys(temp(1)).high_n(1:2,4)));...
    max(abs(Ephys(temp(2)).high_n(1:2,1))) max(abs(Ephys(temp(2)).high_n(1:2,4)));...
  max(abs(Ephys(temp(3)).high_n(1:2,1))) max(abs(Ephys(temp(3)).high_n(1:2,4)));...
   max(abs(Ephys(temp(4)).high_n(1:2,1))) max(abs(Ephys(temp(4)).high_n(1:2,3)))];

cl={'b','k'};
data=[];data=ttx_ipsc;
paired_plot_box(data,cl);
xticklabels({'before','TTX + 4AP'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
xtickangle(45);yticks([0:125:250]); set(gca,'FontSize',10);

cl={'r','k'};
data=[];data=ttx_epsc;
paired_plot_box(data,cl);
xticklabels({'before','TTX + 4AP'});ylabel('EPSC amplitude (pA)');set(gca,'FontSize',10);
xtickangle(45); set(gca,'FontSize',10);
%% TTX Modulation index (not used at moment)
par=[(ttx_ipsc(:,2)-ttx_ipsc(:,1))./(ttx_ipsc(:,2)+ttx_ipsc(:,1)); (ttx_epsc(:,2)-ttx_epsc(:,1))./(ttx_epsc(:,2)+ttx_epsc(:,1))]
s1=[1:4];s2=[5:8]
[statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;
xticklabels({'IPSC' ,'EPSC'});ylabel('TTX modulation index');set(gca,'FontSize',10);xtickangle(45);
%% Show example IPSC TTX before after 
temp=[];temp=find(pyr_washin==1);
cnr=3
ov_min=-20;ov_max=300;
start=4000;
endp=8000;
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
subplot(1,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(start:endp,3),'Color','k','LineWidth',1.5);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_highf(start:endp,2),'Color','b','LineWidth',1.5);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
text(1500,-150,'TTX + 4AP');hold on;text(1500,285,'IPSC before','Color','b');
axis off; set(gca,'FontSize',10);
% subplot(1,2,2)
% plot(Ephys(temp(cnr)).sub_traces_high(start:endp,4),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_highf(start:endp,1),'Color','r','LineWidth',1.2);set(gca,'box','off');
% hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
 %% alternative; just show middle frequency one with median 
 edges = [0:0.25:2];
  fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
hold on;h2=histogram(e_i_ratio_pyr_hf,edges,'Normalization','probability');h2.FaceColor='w';h2.EdgeColor='k';%h2.FaceAlpha=0.5;
 box off;xlabel({'E / I ratio' ; '(25 Hz stim freq.)'});ylabel('Relative counts');
 hold on;plot([1 1],[0 0.4],'--k');
 hold on;plot([nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf))) nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf)))],[0.4 0.4],...
     'Marker','v','MarkerFaceColor','k','MarkerEdgeColor','k');
 %title('25 Hz stim freq.','FontWeight','Normal');
 set(gca,'FontSize',10);
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
xticks([1:1:3]);ylabel('E / I ratio');xticklabels({'1 Hz','25 Hz','50 Hz'});xtickangle(45);set(gca,'FontSize',10);
%% 

%% SPIKING or not PYR vs nonFS and FS
%read out max epsps amplitudes 
epsp_pyr=[];epsp_in=[];
[epsp_pyr] = readout_amp_epsp(Ephys,pyr_k ,2,sr);
[epsp_in] = readout_amp_epsp(Ephys,in_k ,2,sr);

temp=[];maxsF=[];
temp=find(in_k==1);
pv_label=[Ephys(temp).label]==3;
for i=1:length(temp)
    try
    maxsF(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_in{:,i}=Ephys(temp(i)).IV.spikecount;
    stimvec_in{:,i}=Ephys(temp(i)).IV.stimvec;
    catch
    maxsF(i)=NaN;
    spikecount_in{:,i}=NaN;
    stimvec_in{:,i}=NaN;
    end
end
temp=[];maxsF_pyr=[];
temp=find(pyr_k==1);
for i=1:length(temp)
    try
    maxsF_pyr(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_pyr{:,i}=Ephys(temp(i)).IV.spikecount;
    stimvec_pyr{:,i}=Ephys(temp(i)).IV.stimvec;
    catch
    maxsF_pyr(i)=NaN;
    spikecount_pyr{:,i}=NaN;
    stimvec_pyr{:,i}=NaN;
    end
end
%bug in first cell
for t=1:length(Ephys(43).IV.stimvec)
 spike_event_43{:,t}=spike_times(Ephys(43).IV.traces(:,t),1.1);
 spike_log_43(:,t)=~isempty(spike_event_43{:,t});
 spikecount_43(:,t)=length(spike_event_43{:,t});
end
spikecount_in{1, 1}=spikecount_43;
maxsF(1)=max(spikecount_43);
%% Plot histograms for spike frequency 
  edges = [0:8:100];
  fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
  hold on;h2=histogram(maxsF_pyr,edges);h2.FaceColor='k';h2.EdgeColor='k';%h2.FaceAlpha=0.5;
hold on;h2=histogram(maxsF,edges);h2.FaceColor=[0.8 0.8 0.8];h2.EdgeColor='k';%h2.FaceAlpha=0.5;
 box off;xlabel('Max. spike frequency (Hz)');ylabel('Counts');xlim([0 100]);set(gca,'FontSize',10);
legend({'PN','IN'});legend boxoff
%% Classify IN cells based on Fast spiking or not 
temp=[];temp=maxsF;
temp(find(isnan(temp)))=[];
[idx_input_ward, clustering_input, leafOrder] = hca([temp'],0,'ward',2,maxsF',3,0.6);
%% example EPSP traces for PN, nFSIN, FS IN
cnr=16;%
ov_min=-5;ov_max=100;
range=4000:7000;
temp=[];temp=find(pyr_k==1);
fig4=figure;set(fig4, 'Position', [200, 200, 130, 300]);set(gcf,'color','w');
subplot(3,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('PN','Color','k');
axis off;

cnr=3;%
ov_min=-5;ov_max=100;
temp=[];temp=find(in_k==1);
subplot(3,1,2)
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','#A2142F','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('nFS IN','Color','#A2142F');
axis off;

cnr=8;%
ov_min=-5;ov_max=100;
temp=[];temp=find(in_k==1);
subplot(3,1,3)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('FS IN','Color',[0.8500 0.3250 0.0980]);
axis off;
%% showing EPSPs across three gropus
%using a cutoff of 50 Hz
freq_cutoff=50;
a=[];a=find(maxsF>freq_cutoff | pv_label==1);
b=[];b=find(maxsF<freq_cutoff);
c=[];c=find(maxsF<freq_cutoff | pv_label==1);
data=[];data=[epsp_pyr epsp_in(b) epsp_in(a)]';
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 200, 300]);set(gcf,'color','w');
 hold on;scatter(ones(length(epsp_pyr),1),epsp_pyr,25,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
 hold on;scatter(ones(length(epsp_in(b)),1)*2,epsp_in(b),25,'o','MarkerEdgeColor','k','MarkerFaceColor','#A2142F');
 hold on;scatter(ones(length(epsp_in(a)),1)*3,epsp_in(a),25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980]);
%error bars and mean of subthreshold
 temp=[];temp=epsp_in(a);
e=errorbar(1,nanmean(epsp_pyr),nanstd(epsp_pyr)/sqrt(sum(~isnan(epsp_pyr))));
e.Color = 'k';e.CapSize = 10;
  hold on;plot([0.6 1.4],[nanmean(epsp_pyr) nanmean(epsp_pyr)],'Color','k');
e=errorbar(2,nanmean(epsp_in(b)),nanstd(epsp_in(b))/sqrt(sum(~isnan(epsp_in(b)))));
e.Color = '#A2142F';e.CapSize = 10;
  hold on;plot([1.6 2.4],[nanmean(epsp_in(b)) nanmean(epsp_in(b))],'Color','#A2142F');
  e=errorbar(3,nanmean(temp(find(temp<50))),nanstd(temp(find(temp<50)))/sqrt(sum(~isnan(temp(find(temp<50))))));
e.Color = [0.8500 0.3250 0.0980];e.CapSize = 10;
 hold on;plot([2.6 3.4],[nanmean(temp(find(temp<50))) nanmean(temp(find(temp<50)))],'Color',[0.8500 0.3250 0.0980]);
  hold on;p=plot([0 4],[50 50],':k');
 xlim([0 4]);xticks([1:1:3]);ylabel('EPSP amplitude (mV)');xticklabels({'PN','nFS IN','FS IN'});xtickangle(45);
 set(gca,'FontSize',10);
 %ylim([0 35])
 breakyaxis([35 65]);
%% fraction of spikers for PN, nFSIN, FS IN
spike_FS=sum(epsp_in>freq_cutoff)/(sum(maxsF>freq_cutoff));
spike_nFS=0;
spike_PYR=nansum(epsp_pyr>freq_cutoff)/length((epsp_pyr));
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 150, 250]);set(gcf,'color','w');
gr_m=[spike_PYR spike_nFS spike_FS ];  
hold on;
b=bar(1,gr_m(1));b.FaceColor='k';b.FaceAlpha=0.75;
b=bar(2,gr_m(2));b.FaceColor=[0.8500 0.3250 0.0980];b.FaceAlpha=0.25;
b=bar(3,gr_m(3));b.FaceColor=[0.8500 0.3250 0.0980];b.FaceAlpha=1;b.EdgeColor=[0.8500 0.3250 0.0980];
% text(0.7,0.1,num2str(length(e_i_ratio_pyr_long(a1))),'Color','w');hold on;
% text(1.7,0.1,num2str(length(e_i_ratio_pyr_hf(a2))),'Color','w');hold on;
% text(2.7,0.1,num2str(length(e_i_ratio_pyr_hf2(a3))),'Color','w');hold on;
xticks([1:1:3]);ylabel('Fraction of cells firing APs');xticklabels({'PN','nFS IN','FS IN'});xtickangle(45);set(gca,'FontSize',10);
ylim([0 1]);
% % scatter plot
% fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 230, 230]);set(gcf,'color','w');
% scatter(maxsF_pyr,epsp_pyr,25,'o','MarkerEdgeColor','k','MarkerFaceColor','k');hold on;
% scatter(maxsF,epsp_in,25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980]);
% ylabel('EPSP amplitude (mV)');xlabel('max Spike frequency (Hz)');legend({'PYR','IN'}); set(gca,'FontSize',10);legend boxoff
%% compare input to PV and PN: Is input the same or different (under tow conditions: NO TTX 4ap or with TTX 4AP)
%Readout peaks for 
[epsc_pyr ipsc_pyr e_i_ratio_pyr] = readout_amp(Ephys,pyr_cs_pv ,2,2);
[epsc_pv ipsc_pv e_i_ratio_pv] = readout_amp(Ephys,pv_cs_pv ,2,2);
[epsc_pyr_k] = readout_amp(Ephys,pyr_k_pv ,2,1);
[epsc_pv_k] = readout_amp(Ephys,pv_k_pv ,2,1);
%EPSC
data=[];data=[[epsc_pyr epsc_pyr_k]' [epsc_pv epsc_pv_k]'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PN','PV'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);
%IPSC
% data=[];data=[[ipsc_pyr]' [ipsc_pv]'];
% cl={'k',[0.8500 0.3250 0.0980]};
% paired_plot_box(data,cl);hold on;xticklabels({'PN','PV'});ylabel('IPSC Amplitude (pA)');set(gca,'FontSize',10);
% %E/I ratio
% data=[];data=[[e_i_ratio_pyr]' [e_i_ratio_pv]'];
% cl={'k',[0.8500 0.3250 0.0980]};
% paired_plot_box(data,cl);hold on;xticklabels({'PYR','PV'});ylabel('E / I ratio');set(gca,'FontSize',10);
%% Same for cells with TTX
[epsc_pyr_ttx ipsc_pyr_ttx e_i_ratio_pyr_ttx] = readout_amp(Ephys,pyr_ttx ,2,2,1,2);
[epsc_pv_ttx ipsc_pv_ttx e_i_ratio_pv_ttx] = readout_amp(Ephys,pv_ttx ,2,1,2,2);
[epsc_pyr_ttx_sub ipsc_pyr_ttx_sub e_i_ratio_pyr_ttx_sub] = readout_amp(Ephys,pyr_cs_pv_washin,2,2,3,4);
[epsc_pv_ttx_sub ipsc_pv_ttx_sub e_i_ratio_pv_ttx_sub] = readout_amp(Ephys,pv_cs_pv_washin,2,2,3,4);
%replace these accordingly
mm=[];mm=intersect(find(pyr_cs_pv_washin==1),find(pyr_ttx==1));
for i=1:length(mm)
    epsc_pyr_ttx(find(find(pyr_ttx==1)==mm(i)))=epsc_pyr_ttx_sub(find(find(pyr_cs_pv_washin==1)==mm(i)))
end

mm=[];mm=intersect(find(pv_cs_pv_washin==1),find(pv_ttx==1));
for i=1:length(mm)
epsc_pyr_ttx(find(find(pv_ttx==1)==mm(i)))=epsc_pv_ttx_sub(find(find(pv_cs_pv_washin==1)==mm(i)));
end
%EPSC
data=[];data=[epsc_pyr_ttx'  epsc_pv_ttx'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PN','PV'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);


%% 

%% SUPPLEMENTARY %% Cre on vs cre off example cell traces
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
%% Paired comparison EPSC IPSC using the first high frequency pulse for Cre on and Cre off Supp figure b) 
[epsc_on_train ipsc_on_train e_i_ratio_on_train] = readout_amp(Ephys,cre_on_cs ,2,2);
[epsc_off_train ipsc_off_train e_i_ratio_off_train] = readout_amp(Ephys,cre_off_cs ,2,2);
[epsc_on_traink tr trtt] = readout_amp(Ephys,cre_on_k ,2,1);
[epsc_off_traink trr trtttt] = readout_amp(Ephys,cre_off_k ,2,1);
% only using pairs
cl={'m','k'};
data=[];data=[[epsc_on_train epsc_on_traink]' [epsc_off_train epsc_off_traink]'];
paired_plot_box(data,cl);ylabel('Light evoked EPSC (pA)');xticklabels({'Cre+', 'Cre-'});set(gca,'FontSize',10);
data=[];data=[[ipsc_on_train]' [ipsc_off_train]'];
cl={'m','k'};
paired_plot_box(data,cl);ylabel('Light evoked IPSC (pA)');xticklabels({'Cre+', 'Cre-'});set(gca,'FontSize',10);
%% Comparison onset latency for Cre on and Cre off
%CRE ON
temp=[];t_ex_on=[];t_ex_ontr=[];peak_present_on=[];
%temp=find(cre_on_cs==1 | cre_on_k==1);
temp=find(pyr_cs_creon==1);
fc_1=3.5;
close all
for i=1:length(temp)
    try
    if max(Ephys(temp(i)).sub_traces_high(1:1*sr,2))>max(Ephys(temp(i)).sub_traces_high(1:1*sr,1))==1
[t_ex_on(i) t_ex_ontr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,1),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_on(i)=Ephys(temp(i)).high_n(1,1)<0;
    else
      [t_ex_on(i)  t_ex_ontr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,2),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_on(i)=Ephys(temp(i)).high_n(1,2)<0;
    end
    catch
     t_ex_on(i)=NaN;   
    end
end
close all;
%CRE OFF
temp=[];t_ex_off=[]; t_ex_offtr=[];peak_present_off=[];
%temp=find(cre_off_cs==1 | cre_off_k==1);
temp=find(pyr_cs_creoff==1);
fc_1=3.5;
close all
for i=1:length(temp)
    try
    if max(Ephys(temp(i)).sub_traces_high(1:1*sr,2))>max(Ephys(temp(i)).sub_traces_high(1:1*sr,1))==1
    [t_ex_off(i) t_ex_offtr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,1),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_off(i)=Ephys(temp(i)).high_n(1,1)<0;
    else
      [t_ex_off(i) t_ex_offtr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,2),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_off(i)=Ephys(temp(i)).high_n(1,2)<0;
    end
    catch
     t_ex_off(i)=NaN;   
    end
end
close all;
%% using also K-internal cells 
%CRE ON
temp=[];t_ex_kon=[];t_ex_kontr=[];peak_present_kon=[];
%temp=find(cre_on_cs==1 | cre_on_k==1);
temp=find(pyr_k_creon==1);
fc_1=3.5;
close all
for i=1:length(temp)
    try
    if min(Ephys(temp(i)).sub_traces_high(1:1*sr,2))>min(Ephys(temp(i)).sub_traces_high(1:1*sr,1))==1
[t_ex_kon(i) t_ex_kontr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,1),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_kon(i)=Ephys(temp(i)).high_n(1,1)<0;
    else
      [t_ex_kon(i)  t_ex_kontr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,2),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_kon(i)=Ephys(temp(i)).high_n(1,2)<0;
    end
    catch
     t_ex_kon(i)=NaN;   
    end
end
close all;

%CRE OFF
temp=[];t_ex_koff=[]; t_ex_kofftr=[];peak_present_koff=[];
%temp=find(cre_off_cs==1 | cre_off_k==1);
temp=find(pyr_k_creoff==1);
fc_1=3.5;
close all
for i=1:length(temp)
    try
    if min(Ephys(temp(i)).sub_traces_high(1:1*sr,2))>min(Ephys(temp(i)).sub_traces_high(1:1*sr,1))==1
    [t_ex_koff(i) t_ex_kofftr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,1),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_koff(i)=Ephys(temp(i)).high_n(1,1)<0;
    else
      [t_ex_koff(i) t_ex_kofftr(:,i)]=time_to_peak(Ephys(temp(i)).sub_traces_high(:,2),[3000:1:5000],[5000:1:5500],20,0,fc_1);
    peak_present_koff(i)=Ephys(temp(i)).high_n(1,2)<0;
    end
    catch
     t_ex_koff(i)=NaN;   
    end
end
close all;
%% Figure showing onset latencies
cre_on_latency=t_ex_on(find(peak_present_on==1));
cre_on_latency(find(isnan(cre_on_latency)))=[];
cre_off_latency=t_ex_off(find(peak_present_off==1));
cre_off_latency(find(isnan(cre_off_latency)))=[];
cre_kon_latency=t_ex_kon(find(peak_present_kon==1));
cre_kon_latency(find(isnan(cre_kon_latency)))=[];
cre_koff_latency=t_ex_koff(find(peak_present_koff==1));
cre_koff_latency(find(isnan(cre_koff_latency)))=[];
cre_on_latency=[cre_on_latency cre_kon_latency];
cre_off_latency=[cre_off_latency cre_koff_latency];

fig6= figure;set(fig6, 'Name', 'compare latency cre on cre off');set(fig6, 'Position', [200, 300, 150, 250]);set(gcf,'color','w');
gr_m=[nanmean(cre_on_latency) nanmean(cre_off_latency)];  
gr_sem=[nanstd(cre_on_latency)/sqrt(length(cre_on_latency))...
    nanstd(cre_off_latency)/sqrt(length(cre_off_latency))];
hold on;
b=bar(1,gr_m(1));b.FaceColor='k';b.FaceAlpha=0.75;
b=bar(2,gr_m(2));b.FaceColor='k';b.FaceAlpha=0.5;
% hold on;scatter(ones(length(e_i_ratio_pyr_long(a1)),1),e_i_ratio_pyr_long(a1),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf(a2)),1)*2,e_i_ratio_pyr_hf(a2),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf2(a3)),1)*3,e_i_ratio_pyr_hf2(a3),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
hold on;er=errorbar(1:2,gr_m,gr_sem);er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks([1:1:2]);ylabel('Onset latency (ms)');xticklabels({'Cre+','Cre-'});xtickangle(45);set(gca,'FontSize',10);



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


%% Injected current vs Voltage 
a=[];a=find(maxsF>50 | pv_label==1);
b=[];b=find(maxsF<50);
for k=1:25
    for i=1:length(spikecount_in(a))
        try
        s1(i)=spikecount_in{:,a(i)}(k);       
        catch 
         s1(i)=NaN;   
        end
    end
    mean_fs(k)=nanmean(s1);
end
s1=[];
for k=1:25
    for i=1:length(spikecount_in(b))
        try
        s1(i)=spikecount_in{:,b(i)}(k);       
        catch 
         s1(i)=NaN;   
        end
    end
    mean_nfs(k)=nanmean(s1);
end
s1=[];
for k=1:25
    for i=1:length(spikecount_pyr)
        try
        s1(i)=spikecount_pyr{:,i}(k);       
        catch 
         s1(i)=NaN;   
        end
    end
    mean_pyr(k)=nanmean(s1);
end
%% 

 fig4=figure;set(fig4, 'Position', [200, 600, 600, 300]);set(gcf,'color','w');
 for i=1:length(spikecount_in(a))
     p1=plot(stimvec_in{:,a(i)},spikecount_in{:,a(i)},'r', 'MarkerFaceColor','r');p1.Color(4)=3/8;
     hold on;set(gca,'box','off');
 end
hold on;
 for i=1:length(spikecount_in(b))
     p1=plot(stimvec_in{:,b(i)},spikecount_in{:,b(i)},'m', 'MarkerFaceColor','m');p1.Color(4)=3/8;
     hold on;set(gca,'box','off');
 end
hold on;
for i=1:length(spikecount_pyr)
     p1=plot(stimvec_pyr{:,i},spikecount_pyr{:,i},'k', 'MarkerFaceColor','k');p1.Color(4)=3/8;
     hold on;set(gca,'box','off');
end
ylabel('Spike frequency (Hz)');xlabel('Injected current (pA)')
hold on;plot(stimvec_pyr{:,end},mean_fs,'r-o','LineWidth',3);
hold on;plot(stimvec_pyr{:,end},mean_nfs,'m-o','LineWidth',3);
hold on;plot(stimvec_pyr{:,end},mean_pyr,'k-o','LineWidth',3);
%% 
for i=1:length(spikecount_pyr);
pyramidal_cells{:,i}=[stimvec_pyr{:,i}' spikecount_pyr{:,i}'];
end
for i=1:length(spikecount_in(b));
nonfast_interneuron{:,i}=[stimvec_in{:,b(i)}' spikecount_in{:,b(i)}'];
end
for i=1:length(spikecount_in(a));
fast_interneuron{:,i}=[stimvec_in{:,a(i)}' spikecount_in{:,a(i)}'];
end
%% 
data_gain.pyramidal_cells=pyramidal_cells;
data_gain.nonfast_interneuron=nonfast_interneuron;
data_gain.fast_interneuron=fast_interneuron;
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
