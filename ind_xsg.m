function [out] = ind_xsg(pathName)

fac=7;
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
std_bs=fac*std(baseline);
resp_m=max(abs(bs_traces(1000:3000,:)))>std_bs;
for t=1:length(resp_m)
    if resp_m(t)==1
    psc(:,t)=max(abs(bs_traces(1000:3000,t)));
    else resp_m(t)==0;
        psc(:,t)=0;
    end
end
end
f=figure;plot(bs_traces(1:5000,:),'Color',[0.5 0.5 0.5]);
hold on;plot(nanmean(bs_traces(1:5000,:),2),'k');
hold on;plot([base_end*srF (base_end+10)*srF],[50 50],'Color','b','LineWidth',3);%axis off
set(f, 'Position', [200, 200, 250, 230]);ylabel('Synaptic input (pA)')
out=psc';
end