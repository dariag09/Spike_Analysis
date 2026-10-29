% ***********************************************************************************************
% ***********************************************************************************************
% SCRIPT FOR PATCH CLAMP SPIKE ANALYSIS 
%  -- by Daria Garkun
% ***********************************************************************************************
% ***********************************************************************************************

clear all;
close all;

% CONSTANT

AP_Criteria=0;          % threshold parameter to determine an event as an action potential (0 mV)
AP_Threshold=20;        % threshold for AP analysis i.e. threshold for derivative (mV/ms)
DataChannel=3;          % DATA CHANNEL that will be analyzed  

AP_TimeToPeack=1.5;     % window to search AP peak (in milliseconds)
window_before_sweep=4;
window_after_sweep=4;
window_slope_search=20; % window to search rise, decay slopes (in milliseconds)

% DATA LOADING 

[dataFileName,dataPathName] = uigetfile({'*.abf', 'Axon Binary File (*.abf)'},'Select the pClamp file for analysis');
cd(dataPathName);
disp(['Analyzing: ' dataPathName, dataFileName]);                   % print file name and the path to it
[data,si,h] = abfload(dataFileName);                                % loading data from .abf file

% CALCULATING ADDITIONAL VARIABLES FOR FUTURE USE 

sweepsN = size(data,3);                                             % number of sweeps
sweepsDuration=h.dataPtsPerChan*h.si/1000000;                       % duration of sweeps (in seconds)
timeVect=reshape(h.si/1000000:(h.si/1000000):sweepsDuration,[],1);  % time vector
ap_ThresholdToPeack_Window=AP_TimeToPeack*1000/si;                  % window to search for AP peak (in indexes) 

idx_before_sweep=window_before_sweep*1000/si;
idx_after_sweep=window_after_sweep*1000/si;
idx_window_slope_search=window_slope_search*1000/si;

% MEASURE INTERSPIKE INTERVAL 

spikeN=[];

interspikeInterval=NaN(sweepsN,0);
for i = 1:sweepsN 
    Sweep=data(:,DataChannel,i); % reading every sweep separately
    event_Idx = findspike(Sweep, AP_Criteria); % (Sweep, SpkThreshold);
      
    if size(event_Idx,2)>1
        interspikeInterval(i)=min(diff(event_Idx)); 
    else
        interspikeInterval(i)=NaN;            
    end
    spikeN=max([spikeN, size(event_Idx,2)]);    % max spikes
 
end

interspikeInterval_Min=min(interspikeInterval(find(~isnan(interspikeInterval)))); % minimal non-zero and non-NaN element
idx_slope_search=min(idx_window_slope_search, interspikeInterval_Min);

for i = 1:sweepsN 
    
    AP_Threshold_Idx=[]; AP_Threshold_Time=[]; AP_Threshold_Amp=[]; AP_Peak_Idx=[]; AP_Peak_Time=[]; AP_Peak_Amp=[]; clear ('peak');
    AP_rise_slope=[]; AP_decay_slope=[]; AP_Ahp_Voltage=[]; AP_Ahp_Amp=[]; AP_Halfamp_Duration=[];  
   
    Sweep=data(:,DataChannel,i);
    derivative_Sweep = [diff(Sweep);0]*(0.001/(h.si*0.000001));     % (0.001/(h.si*0.000001)) is a conversion  of mV/s to V/s
    event_deriv_Idx=findonset(derivative_Sweep, AP_Threshold);    
        
    disp(['Sweep: ' num2str(i), ' Spike# ' num2str(size(event_deriv_Idx,1))]); 

    if size(event_deriv_Idx(~isnan(event_deriv_Idx)),1)>0

        for j = 1:size(event_deriv_Idx,1) 

            [peak.Amp, peak.LocalIdx]=max(Sweep(event_deriv_Idx(j):(event_deriv_Idx(j)+ap_ThresholdToPeack_Window)));
                     
            if peak.Amp>AP_Criteria        
                AP_Threshold_Idx=[AP_Threshold_Idx; event_deriv_Idx(j)]; 
                AP_Threshold_Time=[AP_Threshold_Time; timeVect(AP_Threshold_Idx(end,1))];
                AP_Threshold_Amp=[AP_Threshold_Amp; Sweep(AP_Threshold_Idx(end,1))];
                
                AP_Peak_Idx=[AP_Peak_Idx; event_deriv_Idx(j)+(peak.LocalIdx)-1];
                AP_Peak_Time=[AP_Peak_Time; timeVect(AP_Peak_Idx(end,1))];
                AP_Peak_Voltage(1:size(Sweep(AP_Peak_Idx,1),1),i)=Sweep(AP_Peak_Idx,1);
                AP_Peak_Amp=[AP_Peak_Amp; AP_Peak_Voltage(end,1)-AP_Threshold_Amp(end,1)];
                
                AP_rise_slope=[AP_rise_slope; max(derivative_Sweep(event_deriv_Idx(j):(event_deriv_Idx(j)+idx_slope_search),1))];
                AP_decay_slope=[AP_decay_slope; min(derivative_Sweep(AP_Peak_Idx(end):(AP_Peak_Idx(end)+idx_slope_search),1))];                
                
                AP_Ahp_Voltage=[AP_Ahp_Voltage; min(Sweep(AP_Peak_Idx(end):(AP_Peak_Idx(end)+idx_slope_search),1))];
                AP_Ahp_Amp=[AP_Ahp_Amp; AP_Ahp_Voltage(end,1)-AP_Threshold_Amp(end,1)];

                AP_Halfamp_Voltage=AP_Peak_Amp(end,1)/2+AP_Threshold_Amp(end,1);
                wave_abovethreshold=find(Sweep(event_deriv_Idx(j):(event_deriv_Idx(j)+idx_slope_search)) > AP_Halfamp_Voltage);
                AP_Halfamp_Duration=[AP_Halfamp_Duration; (wave_abovethreshold(end,1)-wave_abovethreshold(1,1))*si/1000]; % HW in ms
                
            end
            
        end
        
        if (size(event_deriv_Idx(~isnan(event_deriv_Idx)),1)>2) % & (peak.Amp>AP_Criteria)
            
            figure;   

                tiledlayout(3,6);
                nexttile([2 3]);        
                plot (timeVect, Sweep);
                hold on;
                xlim([0.4 2.5]); 
                ylim([-100 50]);
                set(gca,'TickDir','out') % draw the tick marks on the outside
                set(gca,'box','off');
                title(['Sweep # ', num2str(i)]);
                ylabel 'V (mV)';
                xlabel 'Time (s)';

                plot(AP_Threshold_Time(1:size (timeVect(AP_Threshold_Idx,1),1),1), AP_Threshold_Amp(1:size(Sweep(AP_Threshold_Idx,1),1),1), 'x', 'Color', 'g');
                plot(AP_Peak_Time(1:size (timeVect(AP_Peak_Idx,1),1),1), AP_Peak_Voltage(1:size(Sweep(AP_Peak_Idx,1),1),i), 'x', 'Color', 'r');
                hold off;


                nexttile([2 3]);
                plot (Sweep, derivative_Sweep);
                title(['Polar plot, sweep:', num2str(i)]);
                set(gca,'TickDir','out')
                ylabel 'dV';
                xlabel 'V (mV)';

                nexttile;
                plot (AP_Threshold_Amp);
                title('AP Threshold');
                ylabel 'V (mV)';
                xlabel 'AP #';

                nexttile;
                plot (AP_Peak_Amp);
                title('AP Amplitude');
                ylabel 'V (mV)';
                xlabel 'AP #';

                nexttile;
                plot (AP_rise_slope);
                title('Rise slope');
                ylabel 'dV';
                xlabel 'AP #';

                nexttile;
                plot (AP_decay_slope);
                title('Decay slope');
                ylabel 'dV';
                xlabel 'AP #';

                nexttile;
                plot (AP_Ahp_Amp);
                title('AHP amplitude');
                ylabel 'V (mV)';
                xlabel 'AP #';

                nexttile;
                plot (AP_Halfamp_Duration);
                ylim([0 0.4]);
                title('Half width');
                ylabel 'Duration (ms)';
                xlabel 'AP #';

        end

    end
    
 end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function AP_idx = findspike(Wave, EventThreshold)
   clear('AP_idx'); clear('set_cross_plusgi');clear('set_cross_minus');
   
   wave_abovethreshold=find(Wave > EventThreshold); % find indices of samples above the threshold 
    
   if isempty(wave_abovethreshold) < 1  % if there is at least one sample above the threshold
    
       wave_go_up(1)=min(wave_abovethreshold);
       wave_go_down(length(wave_abovethreshold))=max(wave_abovethreshold);

       for i=1:length(wave_abovethreshold)-1
           % detects when there is a discontinuous jump in the
           % indices, thus -> shift from one spike to another. 
           if wave_abovethreshold(i+1) > wave_abovethreshold(i)+1     % index of the right end of the spike interval
               wave_go_up(i+1)=i;
               wave_go_down(i)=i;
           end
       end
       % Code extracts indices of the "above threshold" samples
       % Identifying up and down slopes:
       % Here indices of the truncated "above the threshold" trace are
       % converted into sample indices of the original trace:
       set_cross_plusgi=  wave_abovethreshold(find(wave_go_up));   % find(x) returns nonzero arguments.
       set_cross_minusgi=  wave_abovethreshold(find(wave_go_down));   % find(x) returns nonzero arguments.
       set_cross_minusgi(length(set_cross_plusgi))= wave_abovethreshold(end);
    
       nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows

       for i = 1:nspikes
           % Extract the spike interval
           [~,spikeIdx] = max(Wave(set_cross_plusgi(i):set_cross_minusgi(i)));
           spike_Id(i)=spikeIdx;
           AP_idx(i)=spike_Id(i)+set_cross_plusgi(i)-1;
           % spike_Time(i)=TimeVect((spike_Id(i)+set_cross_plusgi(i)-1),1);
           % spike_Amp(i)=Sweep((spike_Id(i)+set_cross_plusgi(i)-1),3);
       end 
   else 
       AP_idx=NaN; 
      
   end

end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function event_onset = findonset(Wave, EventThreshold)
   clear('AP_idx'); clear('set_cross_plusgi');clear('set_cross_minus');
   
   wave_abovethreshold=find(Wave > EventThreshold); % find indices of samples above the threshold: 
    
   if isempty(wave_abovethreshold) < 1  % if there is at least one sample above the threshold:
    
       % first index posgi = 1st sample above the threshold:
       wave_go_up(1)=min(wave_abovethreshold);
    
       for i=1:length(wave_abovethreshold)-1
           % detects when there is a discontinuous jump in the indices, thus -> shift from one spike to another. 
           if wave_abovethreshold(i+1) > wave_abovethreshold(i)+1     % index of the right end of the spike interval
              wave_go_up(i+1)=i;
           end
       end
       % Code extracts indices of the "above threshold" samples
       % Identifying threshold crossing as wave_go_up

       % Here indices of the truncated "above the threshold" trace are
       % converted into sample indices of the original trace:
       event_onset= wave_abovethreshold(find(wave_go_up));   % find(x) returns nonzero arguments.
       
   else
       event_onset=NaN; 
      
   end

end