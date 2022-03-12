function [norm_traces rise_time decay_time] = constant_traces(Ephys, cell_select,range, range_save, ind1, ind2)

%Ephys: ephys strcuture
%cell_selecter: cells extracted from function cell selecter 
%range: trace range
%range_save: display range 
%ind1, ind2: indices of traces to use 
%% 
sr=20000;
thresh=-50;
temp=[];temp=find(cell_select);
%% 

for i=1:length(temp)
 if   size(Ephys(temp(i)).sub_traces_high(range,:),2)>1 
      
if isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,ind1)<thresh))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,ind2))>max(Ephys(temp(i)).sub_traces_high(range,ind1))

norm_traces(:,i)=Ephys(temp(i)).sub_traces_high(range_save,ind1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,ind1)));
elseif isempty(find(Ephys(temp(i)).sub_traces_high(6.24*sr:end,ind2)<thresh))==0 & ...
        max(Ephys(temp(i)).sub_traces_high(range,ind1))>max(Ephys(temp(i)).sub_traces_high(range,ind2))
norm_traces(:,i)=Ephys(temp(i)).sub_traces_high(range_save,ind2)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,ind2)));
else
    norm_traces(:,i)=ones(length(range_save),1)*NaN;
end
 else
     if sum(Ephys(temp(i)).high_n)==0
         norm_traces(:,i)=ones(length(range_save),1)*NaN;
     else
        norm_traces(:,i)=Ephys(temp(i)).sub_traces_high(range_save,ind1)/abs(min(Ephys(temp(i)).sub_traces_high(range_save,ind1)));
     end
 end
end
%% find time constant 20 / 80% for rise time 
%80 percent of peak
val=-0.8;
for i=1:size(norm_traces,2);
peak_idx=[];
peak_idx=find(norm_traces(:,i)==-1);
[idx_80(i), val_closest(i)] = closest(norm_traces(1:peak_idx(1),i),val);
end
%20 percent of peak
val=-0.2;
for i=1:size(norm_traces,2);
peak_idx=[];
peak_idx=find(norm_traces(:,i)==-1);
[idx_20(i), val_closest(i)] = closest(norm_traces(1:peak_idx(1),i),val);
end
rise_time=(idx_80-idx_20)/(sr/1000);
%% %% find time constant 20 / 80% for decay time 
%80 percent of peak
val=-0.8;
for i=1:size(norm_traces,2);
peak_idx=[];
peak_idx=find(norm_traces(:,i)==-1);
[idx_80(i), val_closest(i)] = closest(norm_traces(peak_idx(1):end,i),val);
end
%20 percent of peak
val=-0.2;
for i=1:size(norm_traces,2);
peak_idx=[];
peak_idx=find(norm_traces(:,i)==-1);
[idx_20(i), val_closest(i)] = closest(norm_traces(peak_idx(1):end,i),val);
end
decay_time=(idx_20-idx_80)/(sr/1000);
end