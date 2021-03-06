function [out]=input_whisker_v1(pathName,maps,tit)
fac=7;
%% Read out individual maps
for i=1:length(maps)
list=dir([char(pathName) 'map' maps{i} filesep '*.xsg']);


load([char(pathName) 'map' maps{i} filesep list.name],'-mat');
sr = header.ephys.ephys.sampleRate;%check sample rate
srF = 1/(1000/sr);
samples_per_sweep = header.ephys.ephys.traceLength*sr;
timebase=1/sr:1/sr:samples_per_sweep/sr; %TR2019: timebase
traces=data.ephys.trace_1;%raw ephys trace
%Filter
cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
type        = 'Butter';
traces = lowpassfilt(traces, order, cutoff, sr, type);
%ind_traces=reshape(traces,[length(traces)/256 256]);
ind_traces=reshape(traces,[length(traces)/128 128]);
base_start=1;
base_end=100;
baseline=ind_traces(base_start*srF:base_end*srF,:);
bs_traces=ind_traces-mean(baseline);%subtract baseline
std_bs=fac*std(baseline);
resp_m=max(abs(bs_traces(1000:3000,:)))>std_bs
for t=1:length(resp_m)
    if resp_m(t)==1
    psc(:,t)=max(abs(bs_traces(1000:3000,t)));
    else resp_m(t)==0;
        psc(:,t)=0;
    end
end
 mapPat=header.mapper.mapper.mapPatternArray;
 array_psc=psc;
  for n=1:numel(mapPat)
  newa(:,find(mapPat==n)) = array_psc(:,n);
  end
  ord_arrays(:,i)=newa;
end
%% multiple repetition criterion
clean_arrays=ord_arrays;
clean_arrays(find(sum(ord_arrays>0,2)<2),:)=0
%Plotting
exc_map=reshape(nanmean(clean_arrays,2),16,8);

if sum(exc_map(:))>0
F = figure;
set(gcf,'color','w');
set(F, 'Name', 'Correlation pial depth');
set(F, 'Position', [200, 200, 180, 230]);

sf=10;    
%define the plot type (2 for excitatory)
explot_type = 2;
%get the pial distance
header.mapper.mapper.soma1Coordinates
map_plot_wV1(exc_map,'',explot_type,F,sf,0,1,header.mapper.mapper.soma1Coordinates);
hold on; title(tit);

else
%% empty matrix
F = figure;
set(gcf,'color','w');
set(F, 'Name', 'Correlation pial depth');
set(F, 'Position', [200, 200, 180, 230]);
imagesc(exc_map);colormap(white);
hold on;
set(gca,'YTick',[1.1, 4.1, 7.1, 9.6, 12.6, 15.1],'YTickLabels',{'L1','L2/3','L4','L5','L6','WM'},...
         'TickLength',[0 0],'XTick',[])
    p1=plot(linspace(0,17,18),2.1.*ones(1,18),':','Color',[0.5 0.5 0.5]);p1.LineWidth=0.25;p1.Color(4) = 0.5;
    p1=plot(linspace(0,17,18),6.1.*ones(1,18),':','Color',[0.5 0.5 0.5]);p1.LineWidth=0.25;p1.Color(4) = 0.5;
    p1=plot(linspace(0,17,18),8.1.*ones(1,18),':','Color',[0.5 0.5 0.5]);p1.LineWidth=0.25;p1.Color(4) = 0.5;
    p1=plot(linspace(0,17,18),11.1.*ones(1,18),':','Color',[0.5 0.5 0.5]);p1.LineWidth=0.25;p1.Color(4) = 0.5;
    p1=plot(linspace(0,17,18),14.1.*ones(1,18),':','Color',[0.5 0.5 0.5]);p1.LineWidth=0.25;p1.Color(4) = 0.5
colorbar;
   x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
  
    adj_x = ((16*69/2)-header.mapper.mapper.soma1Coordinates(1))*(x_lim(2)-x_lim(1))/(16*69);
    adj_y =((16*69/2)-header.mapper.mapper.soma1Coordinates(2))*(y_lim(2)-y_lim(1))/(16*69);
    %plot the center
    plot(adj_x,adj_y,'^k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',2);
    hold on; title(tit);
end
out=clean_arrays;
end