%Reading out xsg files for cracm experiments with S1 ChR2 and recording in
%V1 Munich data from 
%% data storage
%pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0003\map01'
clear all;
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0006\'
maps={'01';'02';'03'};
tit='180313SW0002';
[out]=input_whisker_v1(pathName,maps,tit)
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
%% Baseline subtraction

%% 
array=[];
 newArray=[];
flipFlag = 0;
flipx = 1;
array=bs_traces;
traceLength = size(array,1 );
mapPattern=header.mapper.mapper.mapPatternArray;
if flipx ==1
    mapPattern=fliplr(mapPattern)
end
 for n=1:numel(mapPattern)
 newArray(:,find(mapPattern==n)) = array(:,n);
 end
startTime = 1;
stopTime = traceLength;
showStart = 99*srF;
showStop = 700*srF;
array = newArray(startTime:stopTime, :);
[rows,cols] = size(array);
totalTime = (rows-1)/sr; 
xTimeAxis = linspace(0, totalTime, rows)';
traceAxis = ( 1 : cols );
[sizeX, sizeY] = size(mapPattern);
yFactor = 100; % offset, in pA
%xFactor = 500
scaleBarTxt = 'pA';
if flipFlag == 1
yFactor = -yFactor; 
end
offsetVector = yFactor * ( 0 : cols-1 );
offsetArray = repmat(offsetVector, rows, 1);
array = array-offsetArray;      
% set up the figure -------------------------------------------------------------
x = .14; 
y = .14; 
w = .8; 
h = .8; 
F=figure;
subplotRows = 1; subplotCols = sizeY; plotnum = 0;
for N = 1:sizeY
startInd = (N-1)*sizeX + 1;
endInd = N*sizeX;
plotnum = plotnum+1;

    %     hsub(plotnum) = subplot(subplotRows,subplotCols,plotnum);
    
%     pos1 = 0.025 + (N - 1)*(0.96/sizeY);
%      pos2 = 0.02;
%     pos3 = 0.05;
%     pos4 = 0.96;

   pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    pos2 = 0.02;
    pos3 = 0.038;
    pos4 = 0.96

    hsub(N) = axes('Position', [pos1 pos2 pos3 pos4]);

%     set(gca, 'Position', );

    plot(xTimeAxis(showStart:showStop), array(showStart:showStop,startInd:endInd),'Color','k');
    hold on;
%                 y1=get(gca,'ylim');
%                 x1= redpeak_start/(srF*100);
%                 hold on;
%                 p1=plot([x1 x1],y1,'--','Color','r');
%                 p1.Color(4) = 0.3;
%                 hold on;
% %                 y1=get(gca,'ylim');
% %                 x1=redpeak_end/(srF*100);
% %                 hold on;
% %                 p2=plot([x1 x1],y1,'--','Color','r');
% %                 p2.Color(4) = 0.3;
% %                 hold on;
%                 %blue vertical lines
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_start/(srF*100);
%                 hold on;
%                 p3=plot([x1 x1],y1,'--','Color','b');
%                 p3.Color(4) = 0.3;
%                 hold on;
%                 y1=get(gca,'ylim');
%                 x1=bluepeak_end/(srF*100);
%                 hold on;
%                 p4=plot([x1 x1],y1,'--','Color','b');
%                 p4.Color(4) = 0.3;
    
minval = min(mean(array(1:100,startInd:endInd)));
maxval = max(mean(array(1:100,startInd:endInd)));
tweakFactor = abs(maxval - minval)*0.05;
yrange = [minval-tweakFactor maxval+tweakFactor];
set(gca, 'YLim', yrange);
set(gca, 'XLim', [(showStart-200)/sr (showStop+200)/sr]);
xlabel('Seconds');
set(gcf,'color','w');
%  if flipo==0
%      set(gca, 'XDir','reverse')
% end   
                
end
          
set(findobj(gcf, 'Type', 'axes'), 'Visible', 'off');
    
% scalebar lines
    Y = mean(array(:,end))+yFactor/4;
% Y = min(get(gca, 'YLim'));
    hscalebar = line([.1 .2], [Y Y]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');
    hscalebar = line([.1 .1], [Y Y+yFactor/2]);
    set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');

% scalebar text
    ht(1) = text(.12, Y+yFactor/6, '100 ms'); 
    ht(2) = text(.12, Y+yFactor/3, [num2str(yFactor/2) ' ' scaleBarTxt]); 
    set(ht, 'Color', 'k', 'FontSize', 8, 'Tag', 'scaleBarText');
%axis square;









%% %single traces
clear all
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0007\'
list=dir([char(pathName) filesep '*.xsg']);
for i=1:length(list)
  load([char(pathName) filesep list(i).name],'-mat');  
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
base_start=1;
base_end=100;
baseline=traces(base_start*srF:base_end*srF,:);
bs_traces(:,i)=traces-mean(baseline);%subtract baseline

end
figure;plot(bs_traces)