%Reading out xsg files for cracm experiments with S1 ChR2 and recording in
%V1 Munich data from 
%% data storage
%% single traces c1
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0002\'
c1=ind_xsg(pathName);
%pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0003\map01'
%% c2 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0003\'
maps={'01';'02'};
tit='180302SW0003';
[c2]=input_whisker_v1(pathName,maps,tit);
%% c3 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0006\'
maps={'01';'02';'03'};
tit='180302SW0006';
[c3]=input_whisker_v1(pathName,maps,tit);
%% c4 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0007\'
maps={'01';'02';'03'};
tit='180302SW0007';
[c4]=input_whisker_v1(pathName,maps,tit);
%% c5 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180305iviv\SW0002\'
maps={'01';'02'};
tit='180305SW0002';
[c5]=input_whisker_v1(pathName,maps,tit);
%% c6 single
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180305iviv\SW0003\'
c6=ind_xsg(pathName);

%% c7 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180305iviv\SW0004\'
maps={'03';'04';'05'};
tit='180305SW0004';
[c7]=input_whisker_v1(pathName,maps,tit);
%% c8 maps
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180313iviv\SW0002\'
maps={'01';'02'};
tit='180313SW0002';
[c8]=input_whisker_v1(pathName,maps,tit);
%% single traces c9
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180313iviv\SW0003\'
c9=ind_xsg(pathName);
%% single traces c9
pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180313iviv\SW0004\'
c10=ind_xsg(pathName);
%% all together
all_c=[nanmean(nonzeros(c1));
nanmean(nonzeros(c2))
nanmean(nonzeros(c3))
nanmean(nonzeros(c4))
nanmean(nonzeros(c5))
nanmean(nonzeros(c6))
nanmean(nonzeros(c7))
nanmean(nonzeros(c8))
nanmean(nonzeros(c9))
nanmean(nonzeros(c10))];
all_c(find(isnan(all_c)))=0;
f=figure;b=bar(1:10,all_c);b.FaceColor='w';
xlabel('Cell Number');
ylabel('Average Synaptic Input (pA)')
set(f, 'Position', [200, 200, 250, 230])
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

pathName='C:\Users\slice setup\SimonW\Whisker_V1\data\180302iviv\SW0002\'
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