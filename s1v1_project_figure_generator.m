%% Load data structure for S1 V1 using uipickfiles 
str   = 'D:\Postdoc_Margrie\Projects\Whisker\output_structure';
folder_list = uipickfiles('FilterSpec',str);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%folder where to save figures
save_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Tracing\ManuskriptMT\Revision\review_optophysiology';

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

%PAIRED CS solution without TTX 4AP, only wash in 
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

%%  Read out peaks from train for PN cs for train (long), middle frequency (high), highest (highf)
[epsc_pyr_long ipsc_pyr_long e_i_pyr_train] = readout_amp(Ephys,pyr_cs ,1,2,1,2);
[epsc_pyr_hf ipsc_pyr_hf e_i_ratio_pyr_hf] = readout_amp(Ephys,pyr_cs ,2,2,1,2);
[epsc_pyr_hf2 ipsc_pyr_hf2 e_i_ratio_pyr_hf2] = readout_amp(Ephys,pyr_cs ,3,2,1,2);
%another name for train 
[epsc_pyr_long ipsc_pyr_long e_i_ratio_pyr_long] = readout_amp(Ephys,pyr_cs ,1,2,1,2);
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
%% Paired compariosn of onset latency between EPSC and IPSC  panel d
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
data=[];data=t_ein;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);
%statistics 
%test for normality: 
kstest(t_ex_sub');
kstest(t_in_sub');
%pvalue of paired signrank test:  4.8828e-04
[p1]=signrank(t_ein(:,1),t_ein(:,2));
%for reviewer 3
nanmean(t_ein(:,2))
nanstd(t_ein(:,2))/sqrt(length(t_ein(:,2)))

%% %% TTX modulation index for PN ; not happy with readout; double check 
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
ylim([0 300]);yticks([0:100:300]);

cl={'r','k'};
data=[];data=ttx_epsc;
paired_plot_box(data,cl);
xticklabels({'before','TTX + 4AP'});ylabel('EPSC amplitude (pA)');set(gca,'FontSize',10);
xtickangle(45); set(gca,'FontSize',10);
yticks([0:100:300]);

%% Show example IPSC TTX before after 
temp=[];temp=find(pyr_washin==1);
%before reviews
cnr=3;
ov_min=-20;ov_max=300;
start=4000;
endp=8000;
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
subplot(1,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(start:endp,3),'Color','k','LineWidth',1.5);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_high(start:endp,2),'Color','b','LineWidth',1.5);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
text(1500,-150,'TTX + 4AP');hold on;text(1500,285,'IPSC before','Color','b'); set(gca,'FontSize',10);
%% After reviews
temp=[];temp=find(pyr_washin==1);
cnr=1;
ov_min=-20;ov_max=300;
start=4000;
endp=8000;
fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
subplot(1,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(start:endp,6),'Color','k','LineWidth',1.5);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_high(start:endp,2),'Color','b','LineWidth',1.5);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%% 

% text(1300,-40,'TTX + 4AP (0 mV)');hold on;text(1300,285,'IPSC before (0 mV)','Color','b'); set(gca,'FontSize',10);
% line([-30 470], [-30 -30], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
% line([-30 -30], [-30 20], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
% ax1 = gca;                  
% ax1.YAxis.Visible = 'off';ax1.XAxis.Visible = 'off';ax1.LineWidth=1;hold on;xticks([]);
%% save as pdf 
cd(save_folder);saveas(gcf, 'IPSCexampleTTX_cell13.pdf');
%% 
% temp=[];temp=find(pyr_washin==1);
% cnr=1;
% ov_min=-20;ov_max=300;
% start=4000;
% endp=8000;
% fig4=figure;set(fig4, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
% subplot(1,1,1)
% plot(Ephys(temp(cnr)).sub_traces_train(start:endp,4),'Color','k','LineWidth',1.5);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_train(start:endp,3),'Color','r','LineWidth',1.5);set(gca,'box','off');
% hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
% text(1500,-150,'TTX + 4AP');hold on;text(1500,285,'IPSC before','Color','b'); set(gca,'FontSize',10);
% % subplot(1,2,2)
% % plot(Ephys(temp(cnr)).sub_traces_high(start:endp,4),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
% % hold on;plot(Ephys(temp(cnr)).sub_traces_highf(start:endp,1),'Color','r','LineWidth',1.2);set(gca,'box','off');
% % hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');

 %% E/I ration alternative; just show middle frequency one with median 
% sig_elong =[];sig_ehf = []; sig_ehf2=[];
%  sig_elong+sig_ehf+sig_ehf2

 edges = [0:0.25:2];
  fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
hold on;h2=histogram(e_i_ratio_pyr_hf,edges,'Normalization','probability');h2.FaceColor='w';h2.EdgeColor='k';%h2.FaceAlpha=0.5;
 box off;xlabel({'E / I ratio' ; '(5 Hz stim freq.)'});ylabel('Relative counts');ylim([0 0.6]);
 hold on;plot([1 1],[0 0.5],'--k');
 hold on;plot([nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf))) nanmedian(e_i_ratio_pyr_hf(~isinf(e_i_ratio_pyr_hf)))],[0.5 0.5],...
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
xticks([1:1:3]);ylabel('E / I ratio');xticklabels({'1 Hz','5 Hz','10 Hz'});xtickangle(45);set(gca,'FontSize',10);
%% 

%% SPIKING or not PYR vs nonFS and FS
%read out max epsps amplitudes 
epsp_pyr=[];epsp_in=[];
[epsp_pyr] = readout_amp_epsp(Ephys,pyr_k ,2,sr);
[epsp_in] = readout_amp_epsp(Ephys,in_k ,2,sr);

% epsp_pyr_long=[];epsp_in_long=[];
% [epsp_pyr_long] = readout_amp_epsp(Ephys,pyr_k ,1,sr);
% [epsp_in_long] = readout_amp_epsp(Ephys,in_k ,1,sr);
% 
% epsp_pyr_hf2=[];epsp_in_hf2=[];
% [epsp_pyr_hf2] = readout_amp_epsp(Ephys,pyr_k ,3,sr);
% [epsp_in_hf2] = readout_amp_epsp(Ephys,in_k ,3,sr);

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
maxsF(1)=NaN;
maxsF(11)=NaN;
%% Plot histograms for spike frequency 
  edges = [0:8:100];
  fig4=figure;set(fig4, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
  hold on;h2=histogram(maxsF_pyr,edges,'Normalization','probability');h2.FaceColor='k';h2.EdgeColor='k';%h2.FaceAlpha=0.5;
hold on;h2=histogram(maxsF,edges,'Normalization','probability');h2.FaceColor=[0.8 0.8 0.8];h2.EdgeColor='k';%h2.FaceAlpha=0.5;
 box off;xlabel('Max. spike frequency (Hz)');ylabel('Rleative Counts');xlim([0 100]);set(gca,'FontSize',10);
legend({'PN','IN'});legend boxoff
%% Classify IN cells based on Fast spiking or not using HCA
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
ylim([ov_min-5 12]);title('PN','Color','k');
axis off;

cnr=3;%
ov_min=-5;ov_max=100;
temp=[];temp=find(in_k==1);
subplot(3,1,2)
plot(Ephys(temp(cnr)).sub_traces_high(range,1),'Color','#A2142F','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-5 12]);title('nFS IN','Color','#A2142F');
%axis off;

%8 is in the paper, 19 is good example
cnr=8;%
ov_min=-5;ov_max=100;
temp=[];temp=find(in_k==1);
subplot(3,1,3)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('FS IN','Color',[0.8500 0.3250 0.0980]);
%axis off;
%% 
cd(save_folder);saveas(gcf, 'examples_current_clamp_celltypes.pdf');
%% REVIEWER nr 3, spike latencies questions 
temp=[];temp=find(in_k==1);
spike_temp=[8 9 10 15 16 19];trace_temp=[2 1 2 2 1 2];
for i=1:length(spike_temp)
    trace_pvs1(i,:)=Ephys(temp(spike_temp(i))).sub_traces_high(:,trace_temp(i));
end
for i=1:length(spike_temp)
    trace_pvs2(i,:)=Ephys(temp(spike_temp(i))).sub_traces_highf(:,trace_temp(i));
end
for i=1:length(spike_temp)
    trace_pvs3(i,:)=Ephys(temp(spike_temp(i))).sub_traces_train(:,trace_temp(i));
end
stdfac=10;
time_pv=[];pktime_pv=[];
for i=1:size(trace_pvs,1)
    [time_pv1(i) pktime_pv1(i)]=time_of_resp(trace_pvs1(i,:),1:5000,5001:5500,stdfac); 
end
for i=1:size(trace_pvs,1)
    [time_pv2(i) pktime_pv2(i)]=time_of_resp(trace_pvs2(i,:),1:5000,5001:5500,stdfac); 
end
for i=1:size(trace_pvs,1)
    [time_pv3(i) pktime_pv3(i)]=time_of_resp(trace_pvs3(i,:),1:5000,5001:5500,stdfac); 
end
pktime_allstim=[pktime_pv1/sr*1000; pktime_pv2/sr*1000  ;pktime_pv3/sr*1000];
mean_sp_pvtimes=nanmean(min(pktime_allstim));
sem_sp_pvtimes=nanstd(min(pktime_allstim))/sqrt(length(min(pktime_allstim)));
%% Plot spikes with peak times
range=[];range=4800:5400
fig4=figure;set(fig4, 'Position', [200, 200, 230, 300]);set(gcf,'color','w');
plot(Ephys(temp(8)).sub_traces_train(range,2),'Color',[0.1 0.1 0.1],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,1))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.1 0.1 0.1]);
plot(Ephys(temp(9)).sub_traces_high(range,1),'Color',[0.2 0.2 0.2],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,2))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.2 0.2 0.2]);
plot(Ephys(temp(10)).sub_traces_highf(range,2),'Color',[0.4 0.4 0.4],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,3))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.4 0.4 0.4]);
plot(Ephys(temp(15)).sub_traces_high(range,2),'Color',[0.55 0.55 0.55],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,4))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.55 0.55 0.55]);
plot(Ephys(temp(16)).sub_traces_train(range,1),'Color',[0.7 0.7 0.7],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,5))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.7 0.7 0.7]);
plot(Ephys(temp(19)).sub_traces_train(range,2),'Color',[0.8 0.8 0.8],'LineWidth',1.3);hold on;scatter(min(pktime_allstim(:,6))*sr/1000+200,100,30,'v','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.8 0.8 0.8]);
set(gca,'box','off');
line([200, 200], get(gca, 'ylim'), 'color', ([0 191 255]/256), 'linestyle', ':','LineWidth',1);set(gca,'FontSize',12);hold on;
line([500 600], [-20 -20], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
line([-15 -15], [-20 -10], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
ax1 = gca;                  
ax1.YAxis.Visible = 'off';ax1.XAxis.Visible = 'off';ax1.LineWidth=1;hold on;xticks([]);
%% save as pdf 
cd(save_folder);saveas(gcf, 'spike_latencies_FSIN.pdf');
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
 xlim([0 4]);xticks([1:1:3]);ylabel('Light evoked EPSP amplitude (mV)');xticklabels({'PN','nFS IN','FS IN'});xtickangle(45);
 set(gca,'FontSize',10);
 %ylim([0 35])
 breakyaxis([35 65]);
%% fraction of spikers for PN, nFSIN, FS IN
spike_FS=nansum(epsp_in(a)>freq_cutoff)/(length(a));
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
%% compare input to PV and PN: Is input the same or different (under two conditions: NO TTX 4ap or with TTX 4AP)
%Readout peaks for 
[epsc_pyr ipsc_pyr e_i_ratio_pyr] = readout_amp(Ephys,pyr_cs_pv ,2,2);
[epsc_pv ipsc_pv e_i_ratio_pv] = readout_amp(Ephys,pv_cs_pv ,2,2);
[epsc_pyr_k] = readout_amp(Ephys,pyr_k_pv ,2,1);
[epsc_pv_k] = readout_amp(Ephys,pv_k_pv ,2,1);
%EPSC NO TTX 4AP
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
%for TTX wash in cell use second half of recording where TTX/4AP was present 
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

mm=[];mm=intersect(find(pv_cs_pv_washin==1),find(pv_ttx==1));
for i=1:length(mm)
    ipsc_pyr_ttx(find(find(pyr_ttx==1)==mm(i)))=ipsc_pyr_ttx_sub(find(find(pyr_cs_pv_washin==1)==mm(i)))
end
mm=[];mm=intersect(find(pv_cs_pv_washin==1),find(pv_ttx==1));
for i=1:length(mm)
ipsc_pyr_ttx(find(find(pv_ttx==1)==mm(i)))=ipsc_pv_ttx_sub(find(find(pv_cs_pv_washin==1)==mm(i)));
end
%% Paired plot EPSC PN vs PV with TTX 4AP
data=[];data=[epsc_pyr_ttx'  epsc_pv_ttx'];
cl={'k',[0.8500 0.3250 0.0980]};
paired_plot_box(data,cl);hold on;xticklabels({'PN','PV'});set(gca,'FontSize',10);
title('TTX + 4AP');ylim([0 1000]);
%set(gca,'ytick',[]);
ylabel('EPSC Amplitude (pA)')
%% show example trace without TTX (I think that is fair)
%PN
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
ylim([ov_min-10 ov_max]);title('PN','Color','k');

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
axis off;
% 

%% Same input leads to spiking in one cell but not the other: what ate the underlying reasons 
%The easiest way to analyze the rise time of your EPSP traces is determining the "20-80% rise time". 
% To do so, you measure the amplitude of your EPSP from the preceding baseline 
% (e. g., before the stimulation artefact of your evoked EPSP), then you determine the 20 and 80% of that 
% EPSP amplitude and make a linear fitting between these two points. This is a good estimation of the rise time of 
% PSPs and PSCs in general, regardless of the nature of the signal. 
% This "20-80% rise time" may be comparable to the results of more elaborated 
% mathematical methods such as single exponential fitting and first derivative.
%PN cells under TTX 4AP
range=[];range=1:60000;
cell_select=pyr_ttx;
range_save=5000:7000;
ind1=1;ind2=2;
norm_traces_pyr=[]; rise_time_pyr=[];decay_time_pyr=[];
[norm_traces_pyr rise_time_pyr decay_time_pyr] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);
% use ttx sub after again accordingly
range=[];range=1:60000;
cell_select=pyr_cs_pv_washin;
range_save=5000:7000;
ind1=3;ind2=4;
norm_traces_pyr_sub=[]; rise_time_pyr_sub=[];decay_time_pyr_sub=[];
[norm_traces_pyr_sub rise_time_pyr_sub decay_time_pyr_sub] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);
mm=[];mm=intersect(find(pyr_cs_pv_washin==1),find(pyr_ttx==1));
for i=1:length(mm)
    rise_time_pyr(find(find(pyr_ttx==1)==mm(i)))=rise_time_pyr_sub(find(find(pyr_cs_pv_washin==1)==mm(i)));
    decay_time_pyr(find(find(pyr_ttx==1)==mm(i)))=decay_time_pyr_sub(find(find(pyr_cs_pv_washin==1)==mm(i)));
end
%% PV cells under TTX 4AP
range=[];range=1:60000;
cell_select=pv_ttx;
range_save=5000:7000;
ind1=1;ind2=2;
norm_traces_pv=[]; rise_time_pv=[]; decay_time_pv=[];
[norm_traces_pv rise_time_pv decay_time_pv] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);
range=[];range=1:60000;
cell_select=pv_cs_pv_washin;
range_save=5000:7000;
ind1=3;ind2=4;
norm_traces_pv_sub=[]; rise_time_pv_sub=[]; decay_time_pv_sub=[];
[norm_traces_pv_sub rise_time_pv_sub decay_time_pv_sub] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);
mm=[];mm=intersect(find(pv_cs_pv_washin==1),find(pv_ttx==1));
for i=1:length(mm)
    rise_time_pv(find(find(pv_ttx==1)==mm(i)))=rise_time_pv_sub(find(find(pv_cs_pv_washin==1)==mm(i)));
    decay_time_pv(find(find(pv_ttx==1)==mm(i)))=decay_time_pv_sub(find(find(pv_cs_pv_washin==1)==mm(i)));
end
%% rise and decay time under TTX
p1=[];p2=[];p1=rise_time_pyr;p2=rise_time_pv;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
hold on;xticks([1 2]);xticklabels({'PN','PV'});set(gca,'FontSize',10);
title('TTX + 4AP');ylabel('Rise time (ms)')


%decay
p1=[];p2=[];p1=decay_time_pyr;p2=decay_time_pv;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
hold on;xticks([1 2]);xticklabels({'PN','PV'});set(gca,'FontSize',10);
title('TTX + 4AP');ylabel('Decay time (ms)');
%% Same without TTX 4AP
%rise
range=[];range=1:60000;
cell_select=pyr_cs_pv;
%one is very delayed
range_save=5000:20000;
ind1=1;ind2=2;
norm_traces_pyr_noTTX=[]; rise_time_pyr_noTTX=[];decay_time_pyr_noTTX=[];
[norm_traces_pyr_noTTX rise_time_pyr_noTTX decay_time_pyr_noTTX] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);
%remove outlier
rise_time_pyr_noTTX(3)=[];
decay_time_pyr_noTTX(3)=[];
rise_time_pyr_noTTX(4)=[];
decay_time_pyr_noTTX(4)=[];

range=[];range=1:60000;
cell_select=pv_cs_pv;
range_save=5000:7000;
ind1=1;ind2=2;
norm_traces_pv_noTTX=[]; rise_time_pv_noTTX=[];decay_time_pv_noTTX=[];
[norm_traces_pv_noTTX rise_time_pv_noTTX decay_time_pv_noTTX] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2);

%% rise and decay time without TTX
p1=[];p2=[];p1=rise_time_pyr_noTTX;p2=rise_time_pv_noTTX;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);

p1=[];p2=[];p1=decay_time_pyr_noTTX;p2=decay_time_pv_noTTX;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
%% Pooled together 
g1=[];g2=[];
p1=[];p2=[];p1=[rise_time_pyr rise_time_pyr_noTTX];p2=[rise_time_pv rise_time_pv_noTTX];
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
hold on;xticks([1 2]);xticklabels({'PN','PV'});set(gca,'FontSize',10);
%title('TTX + 4AP');
ylabel('Rise time (ms)')

g1=[];g2=[];
p1=[];p2=[];p1=[decay_time_pyr decay_time_pyr_noTTX];p2=[decay_time_pv decay_time_pv_noTTX];
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
hold on;xticks([1 2]);xticklabels({'PN','PV'});set(gca,'FontSize',10);
%title('TTX + 4AP');
ylabel('Decay time (ms)')
%% GAIN input output
%% Injected current vs Voltage 
a=[];a=find(maxsF>50 | pv_label==1);
b=[];b=find(maxsF<50);
%25 is maximum current steps
for k=1:25
    for i=1:length(spikecount_in(a))
        try
        s1(i)=spikecount_in{:,a(i)}(k);  
        catch 
         s1(i)=NaN;  
        end
        
    end
    mean_fs(k)=nanmean(s1);
    sem_fs(k)=nanstd(s1)/sqrt(length(s1));
end

for i=1:length(spikecount_in(a))
        try 
        max_diff_fs(:,i)=nanmax(diff(spikecount_in{:,a(i)}));
        catch 
         max_diff_fs(:,i)=NaN;
        end
        
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
    sem_nfs(k)=nanstd(s1)/sqrt(length(s1));
end

for i=1:length(spikecount_in(b))
    try 
        max_diff_nfs(:,i)=nanmax(diff(spikecount_in{:,b(i)}));
        catch 
      
         max_diff_nfs(:,i)=NaN;
    end
end


s1=[];max_diff_pyr=[];
for k=1:25
    for i=1:length(spikecount_pyr)
        try
        s1(i)=spikecount_pyr{:,i}(k); 
        catch 
         s1(i)=NaN;
        end
    end
    mean_pyr(k)=nanmean(s1);
    sem_pyr(k)=nanstd(s1)/sqrt(length(s1));
end

for i=1:length(spikecount_pyr)
    try  
        max_diff_pyr(:,i)=nanmax(diff(spikecount_pyr{:,i}));
        catch 
      
         max_diff_pyr(:,i)=NaN;
    end
end
%% Compare slope PN vs PV only

g1=[];g2=[];
p1=[];p2=[];p1=[max_diff_pyr(max_diff_pyr>0)];p2=[max_diff_fs];
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_barplot(par,g1,g2,0);
hold on;xticks([1 2]);xticklabels({'PN','FS'});set(gca,'FontSize',10);
%title('TTX + 4AP');
ylabel('Max Slope');
%% Compare slope PN vs nFS, FS only
gr_nr=3;
g1=[];g2=[];g3=[];
p1=[];p2=[];p3=[];
p1=[max_diff_pyr(max_diff_pyr>0)];p2=[max_diff_nfs];p3=[max_diff_fs]
par=[];par=[p1 p2 p3]';
g1=1:length(p1);
g2=length(p1)+1:length(p1)+1+length(p2)-1;
g3=g2(end)+1:length(par);
gr_nr=3;
color_id={[0 0 0],'#A2142F',[0.8500 0.3250 0.0980]};
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 150, 200]);set(gcf,'color','w');
for i=1:gr_nr
gr_m=[nanmean(par(g1)) nanmean(par(g2)) nanmean(par(g3))];  
gr_sem=[nanstd(par(g1))/sqrt(sum(~isnan(par(g1)))) nanstd(par(g2))/sqrt(sum(~isnan(par(g2))))...
 nanstd(par(g3))/sqrt(sum(~isnan(par(g3))))];
hold on;
b2=bar(i,gr_m(i));b2.FaceColor=color_id{i};
hold on;
plot(ones(1,length(par(g1))),par(g1),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
hold on;
plot(ones(1,length(par(g2)))*2,par(g2),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
hold on;
plot(ones(1,length(par(g3)))*3,par(g3),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
end
hold on;
for i=1:gr_nr
hold on;
er=errorbar(i,gr_m(i),gr_sem(i));er.Color = [0 0 0];er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
end
hold on;xticks([1 2 3]);xticklabels({'PN','nFS IN','FS IN'});set(gca,'FontSize',10);
%title('TTX + 4AP');
ylabel('Max Slope');xtickangle(45);
%Statistical test 
tmp=[];
tmp=par(g3);
tmp(find(isnan(par(g3))))=[];
data_slope=[];data_slope=[par(g1) ;par(g2) ;tmp];
group_id=[];group_id=[ones(length(par(g1)),1); ones(length(par(g2)),1)*2; ones(length(tmp),1)*3];
[p,tbl,stats]  = kruskalwallis(data_slope',group_id');
figure;c = multcompare(stats);
%% Show injected current vs spike frequency
freq_cutoff=50;
a=[];a=find(maxsF>freq_cutoff | pv_label==1);
b=[];b=find(maxsF<freq_cutoff);
c=[];c=find(maxsF<freq_cutoff | pv_label==1);
 fig4=figure;set(fig4, 'Position', [200, 600, 300, 300]);set(gcf,'color','w');
 for i=1:length(spikecount_in(a))
     p1=plot(stimvec_in{:,a(i)},spikecount_in{:,a(i)});p1.Color=[0.8500 0.3250 0.0980];p1.Color(4)=3/8;
    
     hold on;set(gca,'box','off');
 end
hold on;
 for i=1:length(spikecount_in(b))
     p1=plot(stimvec_in{:,b(i)},spikecount_in{:,b(i)});p1.Color='#A2142F';p1.Color(4)=3/8;
     hold on;set(gca,'box','off');
 end

hold on;
for i=1:length(spikecount_pyr)
     p1=plot(stimvec_pyr{:,i},spikecount_pyr{:,i},'k', 'MarkerFaceColor','k');p1.Color(4)=3/8;
     hold on;set(gca,'box','off');
end
ylabel('Spike frequency (Hz)');xlabel('Injected current (pA)')
hold on;p1=plot(stimvec_pyr{:,end},mean_fs,'-o','LineWidth',2);p1.Color=[0.8500 0.3250 0.0980]
 hold on; e=errorbar(stimvec_pyr{:,end},mean_fs,sem_fs);
e.Color = [0.8500 0.3250 0.0980];e.CapSize = 10;
hold on;p1=plot(stimvec_pyr{:,end},mean_nfs,'-o','LineWidth',2);p1.Color='#A2142F';
 hold on; e=errorbar(stimvec_pyr{:,end},mean_nfs,sem_nfs);
e.Color = '#A2142F';e.CapSize = 10;
hold on;p1=plot(stimvec_pyr{:,end},mean_pyr,'-o','LineWidth',2);p1.Color='k';
 hold on; e=errorbar(stimvec_pyr{:,end},mean_pyr,sem_pyr);
e.Color = 'k';e.CapSize = 10;
xlim([0 500]);
% legend({'FS', 'nFS','PN'});
% legend box off

%% Intrinsic properties
[rmp_pyr maxsp_pyr rheo_pyr rin_pyr tau_pyr sag_pyr trace_pyr spike_time_pyr] = passive_readout(Ephys,pyr_k);
[rmp_in maxsp_in rheo_in rin_in tau_in sag_in trace_in spike_time_in] = passive_readout(Ephys,in_k);

%% comparison 
a=[];a=find(maxsF>50 | pv_label==1);
b=[];b=find(maxsF<50);
%Input resistance 
p1=rin_pyr;p2=rin_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_rin]=dual_boxplot(par,g1,g2,0);
ylabel('Input resistance (mOhm)');set(gca,'FontSize',10);
xlim([0 3]);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

%Rheobase
p1=rheo_pyr;p2=rheo_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_rheo]=dual_boxplot(par,g1,g2,0);
xlim([0 3]);ylabel('Threshold current (pA)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

%RMP
p1=rmp_pyr;p2=rmp_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_rmp]=dual_boxplot(par,g1,g2,0);
xlim([0 3]);ylabel('RMP (mV)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);

%tau
p1=tau_pyr;p2=tau_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_tau]=dual_boxplot(par,g1,g2,0);
xlim([0 3]);ylabel('Membrane time constant (ms)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);
%% Active parameters
%first or second spike: first here: 
spike_nr=1;
active_pyr=sp_parameters_pandora(trace_pyr,spike_nr);
active_in=sp_parameters_pandora(trace_in,spike_nr);
 active_pyr(5,:)=active_pyr(5,:)*2;
active_in(5,:)=active_in(5,:)*2;
active_in(15,:)=active_in(15,:)/2;
active_pyr(15,:)=active_pyr(15,:)/2;
%% Subplots for active paramaters 
par1=active_pyr([1 2 3 5 6 7 8 15],:);par2=active_in([1 2 3 5 6 7 8 15],a);
color_id={[0 0 0],[0.8500 0.3250 0.0980]};
str={'APV_{min}(mV)','APV_{peak}(mV)','APV_{thresh}(mV)', 'APVslope_{max} (\DeltamV/\Deltams)','APV_{half} (mV)','APV_{amplitude} (mV)',...
    'AHP_{max}(mV)','AP_{half width} (ms)'};

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 500, 400, 500]);set(gcf,'color','w');
for i=1:size(par1,1)
data={};
data={par1(i,:)' par2(i,:)'};
% gr_m1=nanmean(data{:,1}); 
% gr_m2=nanmean(data{:,2}); 
subplot(3,3,i);hold on;
for k=1:length(data)
    p1=[];p2=[];
    p1=scatter(ones(1,length(data{:,1}))',data{:,1},10,'MarkerEdgeColor',color_id{1});    
    p2=scatter(ones(1,length(data{:,2}))'*2,data{:,2},10,'MarkerEdgeColor',color_id{2}); 
end
hold on;
for m=1:2
er=errorbar(m,nanmean(data{:,m}),nanstd(data{:,m})/sqrt(sum(~isnan(data{:,m}))));er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
er.Color = color_id{m};er.LineWidth=1.5;er.LineStyle = 'none'; 
er.CapSize = 10;hold on;
hold on;plot(m,nanmean(data{:,m}),'o','MarkerFaceColor',color_id{m},'MarkerEdgeColor',color_id{m},'MarkerSize',7);
end
xlim([0 3]);xticks([1 2]);
xticklabels({'PN','FS'});xtickangle(45);
ylabel(str{i});
set(gca,'FontSize',10)
[p k]=ranksum(data{:,1},data{:,2});
    statsout(i)=p;
end
%% difference between rmp and threshold for PN and FS-IN
p1=rmp_pyr;
parcom=[];parcom=[abs(rmp_pyr-par1(3,:))' ;abs(rmp_in(a)-par2(3,:))'];
g1=1:length(p1);
g2=length(p1)+1:length(parcom);
[statsout]=dual_boxplot(parcom,g1,g2,0);
xlim([0 3]);ylabel('\DeltaRMP - APV_{thresh}');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS IN'});xtickangle(45);
%% 

p1=active_pyr(5,:);p2=active_in(5,:);
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_boxplot(par,g1,g2,0);
xlim([0 3]);ylabel('APVslope_{max} (\DeltamV/\Deltams)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS'});
xtickangle(45);

%% 

 Property = ["RMP";"R_{in}";"Tau";"Threshold current";"APV_{min}";"APV_{peak}";"APV_{thresh}";"APVslope_{max}";...
     "APV_{half}";"APV_{amplitude}";"AHP_{max}";"AP_{half width}"];
 %active_pyr([1 2 3 5 6 7 8 15],:)
 PN = [nanmean(rmp_pyr);nanmean(rin_pyr);nanmean(tau_pyr);nanmean(rheo_pyr);nanmean(active_pyr(1,:));nanmean(active_pyr(2,:));...
     nanmean(active_pyr(3,:));nanmean(active_pyr(5,:));nanmean(active_pyr(6,:));nanmean(active_pyr(7,:));nanmean(active_pyr(8,:));nanmean(active_pyr(15,:))];

 act_idx=[1 2 3 5 6 7 8 15];
for i=1:length(act_idx)
    temp3=[];temp3=active_pyr(act_idx(i),:);
     sem_pyr_all(i)=nanstd(temp3)/sqrt(sum(~isnan(temp3)));
end

for i=1:length(act_idx)
    temp3=[];temp3=active_in(act_idx(i),a);
     sem_in_all(i)=nanstd(temp3)/sqrt(sum(~isnan(temp3)));
end

 SEM_PN=[nanstd(rmp_pyr)/sqrt(sum(~isnan(rmp_pyr)));nanstd(rin_pyr)/sqrt(sum(~isnan(rin_pyr)))...
     ;nanstd(tau_pyr)/sqrt(sum(~isnan(tau_pyr)));nanstd(rheo_pyr)/sqrt(sum(~isnan(rheo_pyr)));sem_pyr_all'];

 FS_IN = [nanmean(rmp_in(a));nanmean(rin_in(a));nanmean(tau_in(a));nanmean(rheo_in(a));nanmean(active_in(1,a));nanmean(active_in(2,a));...
     nanmean(active_in(3,a));nanmean(active_in(5,a));nanmean(active_in(6,a));nanmean(active_in(7,a));nanmean(active_in(8,a));nanmean(active_in(15,a))];
 
  SEM_IN=[nanstd(rmp_in(a))/sqrt(sum(~isnan(rmp_in(a))));nanstd(rin_in(a))/sqrt(sum(~isnan(rin_in(a))))...
     ;nanstd(tau_in(a))/sqrt(sum(~isnan(tau_in(a))));nanstd(rheo_in(a))/sqrt(sum(~isnan(rheo_in(a))));sem_in_all'];

 
 pValue = [statsout_rmp;statsout_rin;statsout_tau;statsout_rheo;statsout'];

% 

 intrprop = table(Property,PN,SEM_PN,FS_IN,SEM_IN,pValue);




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
cnr=7;%210907SW001
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
[epsc_on_train ipsc_on_train e_i_ratio_on_train] = readout_amp(Ephys,cre_on_cs ,2,2,1,2);
[epsc_off_train ipsc_off_train e_i_ratio_off_train] = readout_amp(Ephys,cre_off_cs ,2,2,1,2);
[epsc_on_traink tr trtt] = readout_amp(Ephys,cre_on_k ,2,1,1,2);
[epsc_off_traink trr trtttt] = readout_amp(Ephys,cre_off_k ,2,1,1,2);
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
b=bar(1,gr_m(1));b.FaceColor='m';b.FaceAlpha=0.75;
b=bar(2,gr_m(2));b.FaceColor='k';b.FaceAlpha=0.75;
% hold on;scatter(ones(length(e_i_ratio_pyr_long(a1)),1),e_i_ratio_pyr_long(a1),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf(a2)),1)*2,e_i_ratio_pyr_hf(a2),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
% hold on;scatter(ones(length(e_i_ratio_pyr_hf2(a3)),1)*3,e_i_ratio_pyr_hf2(a3),7,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
hold on;er=errorbar(1:2,gr_m,gr_sem);er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
xticks([1:1:2]);ylabel('Onset latency (ms)');xticklabels({'Cre+','Cre-'});xtickangle(45);set(gca,'FontSize',10);

%% REVIEW LOW LASER INTENSITIES 
pn=cell_selecter(Ephys,'label',5,'sol',1,'layer',3,'drugs',2);
pv=cell_selecter(Ephys,'label',4,'sol',1,'layer',3,'drugs',2);

temp=[];temp=find(pn==1);
temp2=[];temp2=find(pv==1);
%% 
%low intensity ramp, read out for PN and PVs
%PN
a1=[];a1=nanmean(Ephys(temp(1)).ladder_p(:,4:end),2);
a2=[];a2=nanmean(Ephys(temp(2)).ladder_p(:,5:7),2);
a3=[];a3=nanmean(Ephys(temp(3)).ladder_p(:,1:3),2);
a4=[];a4=nanmean(Ephys(temp(4)).ladder_p(:,6:10),2);
a5=[];a5=nanmean(Ephys(temp(5)).ladder_p(:,1:4),2);
%PV
b1=[];b1=nanmean(Ephys(temp2(1)).ladder_p(:,1:2),2);
b2=[];b2=nanmean(Ephys(temp2(2)).ladder_p(:,6:10),2);
b3=[];b3=nanmean(Ephys(temp2(3)).ladder_p(:,1:3),2);
b4=[];b4=nanmean(Ephys(temp2(4)).ladder_p(:,5:7),2);
b5=[];b5=nanmean(Ephys(temp2(5)).ladder_p(:,2:7),2);
b6=[];b6=nanmean(Ephys(temp2(6)).ladder_p(:,5:8),2);

%combine
pn_meps=[a1 a2 a3 a4 a5];
pv_meps=[b1 b2 b3 b4 b5 b6];




%% 


fig6= figure;set(fig6, 'Name', 'Review minimal epsps');set(fig6, 'Position', [200, 300, 350, 250]);set(gcf,'color','w');
for i=1:size(pn_meps,2)
plot(pn_meps(:,i),'-','Color',[0 0 0 0.3]);
min_resppn(i)=min(nonzeros(pn_meps(:,i)))
hold on

end
hold on;


for i=1:size(pv_meps,2)
plot(pv_meps(:,i),'-','Color',[0.8500 0.3250 0.0980 0.3]);
min_resppv(i)=min(nonzeros(pv_meps(:,i)))
hold on
end

hold on; e1=errorbar(1:11,nanmean(pv_meps,2),nanstd(pv_meps,[],2)/sqrt(length(pv_meps)),'LineWidth',2)
e1.Color = [0.8500 0.3250 0.0980];e.CapSize = 10;
hold on; e2=errorbar(1:11,nanmean(pn_meps,2),nanstd(pn_meps,[],2)/sqrt(length(pn_meps)),'LineWidth',2)
e2.Color = [0 0 0];e.CapSize = 10;
box off;
xlabel('473 nm intensity steps');ylabel('Light evoked EPSP amplitude (mV)');set(gca,'FontSize',10);
xticks([1:1:11]);
hold off
legend([e1 e2],{'PV','PN'});legend boxoff

%% 

p1=min_resppn;p2=min_resppv;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_rin]=dual_boxplot(par,g1,g2,0);
xticklabels({'PN',' PV'});xtickangle(45);ylabel('minimal light evoked EPSP (mV)');set(gca,'FontSize',10);
% fig6= figure;set(fig6, 'Name', 'Review minimal epsps');set(fig6, 'Position', [200, 300, 350, 250]);set(gcf,'color','w');
% for i=1:size(pn_meps,2)
% plot(pn_meps(:,i)/max(pn_meps(:,i)),'-','Color',[0 0 0 0.3]);
% hold on
% 
% end
% hold on;
% 
% 
% for i=1:size(pv_meps,2)
% plot(pv_meps(:,i)/max(pv_meps(:,i)),'-','Color',[0.8500 0.3250 0.0980 0.3]);
% hold on
% end
% 
% hold on; e1=errorbar(1:11,nanmean(pv_meps./max(pv_meps),2),nanstd(pv_meps./max(pv_meps),[],2)/sqrt(length(pv_meps)),'LineWidth',2)
% e1.Color = [0.8500 0.3250 0.0980];e.CapSize = 10;
% hold on; e2=errorbar(1:11,nanmean(pn_meps./max(pn_meps),2),nanstd(pn_meps./max(pn_meps),[],2)/sqrt(length(pn_meps)),'LineWidth',2)
% e2.Color = [0 0 0];e.CapSize = 10;
% box off;
% xlabel('473 nm intensity steps');ylabel('Normalized amplitude');set(gca,'FontSize',10);
% xticks([1:1:11]);
% hold off
% legend([e1 e2],{'PV','PN'});legend boxoff


%% 

hold on;plot(nanmean(Ephys(temp(2)).ladder_p(:,5:7),2),'-or')
hold on;plot(nanmean(Ephys(temp(3)).ladder_p(:,1:3),2),'-or')
hold on;plot(nanmean(Ephys(temp(4)).ladder_p(:,6:10),2),'-or')
hold on;plot(nanmean(Ephys(temp(5)).ladder_p(:,1:4),2),'-or')

hold on;plot(nanmean(Ephys(temp2(1)).ladder_p(:,1:2),2),'-ob')
hold on;plot(nanmean(Ephys(temp2(2)).ladder_p(:,6:10),2),'-ob')
hold on;plot(nanmean(Ephys(temp2(3)).ladder_p(:,1:3),2),'-ob')
hold on;plot(nanmean(Ephys(temp2(4)).ladder_p(:,5:7),2),'-ob')
hold on;plot(nanmean(Ephys(temp2(5)).ladder_p(:,2:7),2),'-ob')
hold on;plot(nanmean(Ephys(temp2(6)).ladder_p(:,5:8),2),'-ob')



% % %high intensity ramp
% figure;
% plot(nanmean(Ephys(temp(1)).ladder_p(:,4:end),2),'r')
% hold on;plot(nanmean(Ephys(temp(2)).ladder_p(:,8:end),2),'r')
% hold on;plot(nanmean(Ephys(temp(3)).ladder_p(:,4:end),2),'r')
% 
% hold on;plot(nanmean(Ephys(temp2(1)).ladder_p(:,3:end),2),'b')
% hold on;plot(nanmean(Ephys(temp2(2)).ladder_p(:,11:end),2),'b')
% hold on;plot(nanmean(Ephys(temp2(3)).ladder_p(:,4:end),2),'b')
% hold on;plot(nanmean(Ephys(temp2(4)).ladder_p(:,8:end),2),'b')




%% Stuff for Vahid
a=[];a=find(maxsF>50 | pv_label==1);
b=[];b=find(maxsF<50);
for i=1:length(spikecount_pyr);
pyramidal_cells{:,i}=[stimvec_pyr{:,i}' spikecount_pyr{:,i}'];
end
for i=1:length(spikecount_in(b));
nonfast_interneuron{:,i}=[stimvec_in{:,b(i)}' spikecount_in{:,b(i)}'];
end
for i=1:length(spikecount_in(a));
fast_interneuron{:,i}=[stimvec_in{:,a(i)}' spikecount_in{:,a(i)}'];
end
data_gain.pyramidal_cells=pyramidal_cells;
data_gain.nonfast_interneuron=nonfast_interneuron;
data_gain.fast_interneuron=fast_interneuron;

%% Read out Rheobase for Vahid
temp=[];temp=find(pyr_k==1);
for i=1:length(temp)
   
    if ~isempty(Ephys(temp(i)).Rheobase)==1
rheo_spikecount_pn(i,:)=Ephys(temp(i)).Rheobase.spikecount(1:20);
rheo_current_pn(i,:)=Ephys(temp(i)).Rheobase.stimvec(1:20);
    else 
        rheo_spikecount_pn(i,:)=ones(1,20)*NaN;
        rheo_current_pn(i,:)=ones(1,20)*NaN;
    end
end
%% 
temp=[];temp=find(pyr_k==1);
for i=1:length(temp)
   
    if ~isempty(Ephys(temp(i)).Rheobase)==1
rheo_spikecount_pn(i,:)=Ephys(temp(i)).Rheobase.spikecount(1:20);
rheo_current_pn(i,:)=Ephys(temp(i)).Rheobase.stimvec(1:20);
    else 
        rheo_spikecount_pn(i,:)=ones(1,20)*NaN;
        rheo_current_pn(i,:)=ones(1,20)*NaN;
    end
end
%% 
temp=[];temp=find(in_k==1);
for i=1:length(temp)
   
    if ~isempty(Ephys(temp(i)).Rheobase)==1
rheo_spikecount_in(i,:)=Ephys(temp(i)).Rheobase.spikecount(1:20);
rheo_current_in(i,:)=Ephys(temp(i)).Rheobase.stimvec(1:20);
    else 
        rheo_spikecount_in(i,:)=ones(1,20)*NaN;
        rheo_current_in(i,:)=ones(1,20)*NaN;
    end
end

%% 
rheobase.pn_current=rheo_current_pn;
rheobase.pn_spikecount=rheo_spikecount_pn;
rheobase.nonfast_current=rheo_current_in(b,:);
rheobase.nonfast_spikecount=rheo_spikecount_in(b,:);
rheobase.fast_current=rheo_current_in(a,:);
rheobase.fast_spikecount=rheo_spikecount_in(a,:);
%% 




% %% Injected current vs Voltage 
% a=[];a=find(maxsF>50 | pv_label==1);
% b=[];b=find(maxsF<50);
% for k=1:25
%     for i=1:length(spikecount_in(a))
%         try
%         s1(i)=spikecount_in{:,a(i)}(k);       
%         catch 
%          s1(i)=NaN;   
%         end
%     end
%     mean_fs(k)=nanmean(s1);
% end
% s1=[];
% for k=1:25
%     for i=1:length(spikecount_in(b))
%         try
%         s1(i)=spikecount_in{:,b(i)}(k);       
%         catch 
%          s1(i)=NaN;   
%         end
%     end
%     mean_nfs(k)=nanmean(s1);
% end
% s1=[];
% for k=1:25
%     for i=1:length(spikecount_pyr)
%         try
%         s1(i)=spikecount_pyr{:,i}(k);       
%         catch 
%          s1(i)=NaN;   
%         end
%     end
%     mean_pyr(k)=nanmean(s1);
% end
% %% 
% 
%  fig4=figure;set(fig4, 'Position', [200, 600, 600, 300]);set(gcf,'color','w');
%  for i=1:length(spikecount_in(a))
%      p1=plot(stimvec_in{:,a(i)},spikecount_in{:,a(i)},'r', 'MarkerFaceColor','r');p1.Color(4)=3/8;
%      hold on;set(gca,'box','off');
%  end
% hold on;
%  for i=1:length(spikecount_in(b))
%      p1=plot(stimvec_in{:,b(i)},spikecount_in{:,b(i)},'m', 'MarkerFaceColor','m');p1.Color(4)=3/8;
%      hold on;set(gca,'box','off');
%  end
% hold on;
% for i=1:length(spikecount_pyr)
%      p1=plot(stimvec_pyr{:,i},spikecount_pyr{:,i},'k', 'MarkerFaceColor','k');p1.Color(4)=3/8;
%      hold on;set(gca,'box','off');
% end
% ylabel('Spike frequency (Hz)');xlabel('Injected current (pA)')
% hold on;plot(stimvec_pyr{:,end},mean_fs,'r-o','LineWidth',3);
% hold on;plot(stimvec_pyr{:,end},mean_nfs,'m-o','LineWidth',3);
% hold on;plot(stimvec_pyr{:,end},mean_pyr,'k-o','LineWidth',3);
%% 

% %% TTX wash in 
% all_cs_ttx = cell_selecter(Ephys,'sol',2,'drugs',1);
% temp=[];
% temp=find(all_cs_ttx==1);
% cnr=3
% ov_min=-20;ov_max=300;
% fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
% subplot(1,2,1)
% plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,3),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,2),'Color','b','LineWidth',1.2);set(gca,'box','off');
% hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
% subplot(1,2,2)
% plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,4),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_highf(1:1*sr,1),'Color','r','LineWidth',1.2);set(gca,'box','off');
% hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
% %% IPSC quantification
% cl={'b',[0.5 0.5 0.5]};
% ttx_ipsc=[max(abs(Ephys(temp(3)).highf_p(:,2))) max(abs(Ephys(temp(3)).highf_p(:,3)));...
%     max(abs(Ephys(temp(2)).highf_p(:,2))) max(abs(Ephys(temp(2)).highf_p(:,3)));...
%   max(abs(Ephys(temp(1)).highf_p(:,1))) max(abs(Ephys(temp(1)).highf_p(:,3))) ];
% data=[];data=ttx_ipsc;
% fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
% hold on;
% for i=1:length(data)
%      pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
% end
% hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
%     'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
% box off;set(gca,'FontSize',10);
%  [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
% title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
% xticklabels({'no TTX','TTX'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
% %% EPSC
% ttx_epsc=[max(abs(Ephys(temp(3)).highf_n(:,1))) max(abs(Ephys(temp(3)).highf_n(:,4)));...
%     max(abs(Ephys(temp(2)).highf_n(:,1))) max(abs(Ephys(temp(2)).highf_n(:,4)));...
%   max(abs(Ephys(temp(1)).highf_n(:,3))) max(abs(Ephys(temp(1)).highf_n(:,4))) ];
% data=[];data=ttx_epsc;
% fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
% hold on;
% for i=1:length(data)
%      pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
% end
% hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
%     'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
% box off;set(gca,'FontSize',10);
%  [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
% title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
% xticklabels({'no TTX','TTX'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
% %% TTX modulation index
% par=[(ttx_ipsc(:,2)-ttx_ipsc(:,1))./(ttx_ipsc(:,2)+ttx_ipsc(:,1)); (ttx_epsc(:,2)-ttx_epsc(:,1))./(ttx_epsc(:,2)+ttx_epsc(:,1))]
% s1=[1:3];s2=[4:6]
% [statsout]=dual_barplot(par,s1,s2,2);xticks([1:1:2]);hold on;xticklabels({'IPSC' ,'EPSC'});ylabel('TTX modulation index');set(gca,'FontSize',10);xtickangle(45);
% %% 
% all_red=[];all_nred=[];
% all_red = cell_selecter(Ephys,'label',1,'pair',1,'geno',7);
% all_nred = cell_selecter(Ephys,'label',0,'pair',1,'geno',7);
% %% 
% all_red=[];all_nred=[];
% all_red = cell_selecter(Ephys,'label',1,'geno',7);
% all_nred = cell_selecter(Ephys,'label',0,'geno',7);
% %% 
% red_epsc=[];red_epscf=[];nred_epsc=[];nred_epscf=[];
% temp=[];
% for i=1:length(find(all_red==1));
%    temp=find(all_red==1);
%    red_epsc(i)=max(abs(Ephys(temp(i)).train_n(:)));
% end
% temp=[];
%    for i=1:length(find(all_nred==1));
%    temp=find(all_nred==1);
%    nred_epsc(i)=max(abs(Ephys(temp(i)).train_n(:)));
%    end
%    
% temp=[];
% for i=1:length(find(all_red==1));
%    temp=find(all_red==1);
%    red_epscf(i)=max(abs(Ephys(temp(i)).high_n(:)));
% end
% temp=[];
%    for i=1:length(find(all_nred==1));
%    temp=find(all_nred==1);
%    nred_epscf(i)=max(abs(Ephys(temp(i)).high_n(:)));
%    end
%    
% %    temp=[];
% % for i=1:length(find(all_red==1));
% %    temp=find(all_red==1);
% %    red_epschf(i)=max(abs(Ephys(temp(i)).highf_n(:)));
% % end
% % temp=[];
% %    for i=1:length(find(all_nred==1));
% %    temp=find(all_nred==1);
% %    nred_epschf(i)=max(abs(Ephys(temp(i)).highf_n(:)));
% %    end
%    
%    %% All against all
%    par=[red_epsc nred_epsc]
%    s1=1:length(red_epsc)
%    s2=length(red_epsc)+1:length(nred_epsc)+length(red_epsc)
%   [statsout]=dual_barplot(par,s1,s2,0);xticks([1:1:2]);hold on;xticklabels({'Cre+' ,'Cre-'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);xtickangle(45);
% 
% %% only using pairs
% cl={'r',[0.5 0.5 0.5]};
% data=[];data=[red_epscf' nred_epscf'];
% fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
% hold on;
% for i=1:length(data)
%      pl=plot([1,2],[data(:,1),data(:,2)],'color','k');    
% end
% 
% hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
%     'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
% %hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
% %hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
% box off;set(gca,'FontSize',10);
% %hold on;errorbar([0.75 2.25],nanmean(data),nanstd(data,[],1)/sqrt(length(data)),'ok','MarkerFaceColor','r','Markersize',7);
% hold on;errorbar([0.75],nanmean(data(:,1)),nanstd(data(:,1),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor','r','Markersize',7);
% hold on;errorbar([2.25],nanmean(data(:,2)),nanstd(data(:,2),[],1)/sqrt(length(data(:,1))),'ok','MarkerFaceColor',[0.5 0.5 0.5],'Markersize',7);
%  [p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
% title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
% xticklabels({'Cre+','Cre-'});ylabel('EPSC Amplitude (pA)');set(gca,'FontSize',10);
% 
% %% 
% 
% cnr=6
% ov_min=-400;ov_max=600;
% temp=[];
% temp=find(all_red==1);
% fig4=figure;set(fig4, 'Position', [200, 200, 600, 300]);set(gcf,'color','w');
% subplot(1,2,1)
% plot(Ephys(temp(cnr)).sub_traces_high(1:2*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
% ylim([-200 20])
% ov_min=-400;ov_max=600;
% title('Cre+')
% temp=[];
% temp=find(all_nred==1);
% subplot(1,2,2)
% 
% plot(Ephys(temp(cnr)).sub_traces_high(1:2*sr,1),'Color',[0.5 0.5 0.5],'LineWidth',1);set(gca,'box','off');
% ylim([-200 20])
% title('Cre-')
% %% 
% 
%  plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
% hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);   
% hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
% % %% 


%% 
temp=[];temp=find(pyr_k==1);
for i=1:length(temp)
    if ~isempty(Ephys(temp(i)).Rheobase)==1
    rheo_fix_pk(i)=Ephys(temp(i)).Rheobase.rheo;
    else ~isempty(Ephys(temp(i)).Rheobase)==0
        if ~isempty(Ephys(temp(i)).IV)==1
        tr=[];tr=find(Ephys(temp(i)).IV.spikecount>0);
        rheo_fix_pk(i)=Ephys(temp(i)).IV.stimvec(tr(1));
        else
        rheo_fix_pk(i)=NaN;
        end
    end
end
%% 

p1=rheo_fix_pk;p2=rheo_in;
par=[];par=[p1 p2(a)]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout]=dual_boxplot(par,g1,g2,0);
xlim([0 3]);ylabel('Threshold current (pA)');set(gca,'FontSize',10);xticks([1 2]);xticklabels({'PN','FS'});xtickangle(45);
%% 

