function z=fft_comparison_average_pulses_envelope_v2(filename)

close all;
    %raw data extraction
[data,Fs]=audioread(filename);
    %gp variables and parameters
time_experiment=180;
max_index=time_experiment*Fs;
pulse_length=(10/1000)*Fs;

threshold=1/1000;     %ideally it should be 1mV but to exclude reflections 2.5mV makes more sense
threshold_back_samples=0;

samples=Fs;
freq_set=((-samples/2):1:(samples/2)-1)*(Fs/samples);

total_number_of_pulses=0;
pulses_with_valid_freq=0;
last_index=(-300*Fs/1000);  %in order to avoid missing first pulse
average_envelope_freq=1;
no_of_sets=50000;
data_extracted=zeros(samples,1,1);
single_pulse=zeros(samples,1);
alternate_colour_flag=1;
%fft_set=zeros(samples,1,no_of_sets);
%fft_temp=zeros(samples,1);
    %raw data plot
mean_value=mean(data);
data=data-mean_value;
figure
time=0:1/Fs:time_experiment;
time=time(1:46080000);
p1=plot(time,data(1:max_index),'b');
xlabel('time(s)');
ylabel('Ampl.');
    %filter and filtered data plot
fc = 69000;
[b,a] = butter(6,fc/(Fs/2));
dataFilt = filter(b,a,data);
hold on
p2=plot(time,dataFilt(1:max_index),'r');
h = [p1(1);p2];
legend(h,'Raw data','Filtered'); 
    %extracting pulses from data and plot
index=1;
    while( index < max_index)
        index=index+1;
        if ( dataFilt(index) > (threshold) )
            %%
            lower_limit=index-threshold_back_samples;
            if(lower_limit<0)
               index=index+3*pulse_length;
            else
                %%
                total_number_of_pulses=total_number_of_pulses+1;
                upper_limit=lower_limit+pulse_length-1;
                if upper_limit > max_index
                    break;
                end
                data_extracted(1:pulse_length,1,1)=dataFilt(lower_limit:upper_limit);
                single_pulse(:,1)=data_extracted(:,1,1);
                [fft_temp,temp_f_peak]=fft_pulse(single_pulse,samples,freq_set);
                if temp_f_peak > 66500
                   % index
                   % index*1/Fs
                   % last_index
                   % last_index*1/Fs
                   % index-last_index
                   % ((index-last_index)/Fs)*1000
                   % (3/1000)
                    %%
                    if((index-last_index)/Fs>((300/1000)))
                        %%
                        pulses_with_valid_freq=pulses_with_valid_freq+1;
                        %fft_set=fft_temp(samples,1,pulses_with_valid_freq);
                        extracted_peak(pulses_with_valid_freq,1)=temp_f_peak;
                        if pulses_with_valid_freq>1 && mod(pulses_with_valid_freq,8)==0
                            extracted_average_peak(average_envelope_freq,1)=mean(extracted_peak(pulses_with_valid_freq-7:pulses_with_valid_freq,1));
                            average_envelope_freq=average_envelope_freq+1;
                        else
                        end
                        if mod(pulses_with_valid_freq,8)==1
                            alternate_colour_flag=0;
                        else
                            alternate_colour_flag=1;
                        end
                            %plot only valid pulses, this will help in next
                            %version of program/algo
                        time_for_one_pulse=0:1/Fs:1;
                        time_for_one_pulse=(lower_limit*(1/Fs))+time_for_one_pulse(1:256000);
                        hold on
                        if (alternate_colour_flag==1)
                            plot(time_for_one_pulse, data_extracted(:,1,1),'g');
                        else
                            plot(time_for_one_pulse, data_extracted(:,1,1),'k');
                        end
                            last_index=index;       %used for rejections of echos....
                            index=index+3*pulse_length;
                    else
                        index=index+3*pulse_length; %in case of echo
                    end
                else
                        index=index+(1/1000)*Fs; %1ms in case of invalid pulse 
                end
                %if total_number_of_pulses >no_of_sets
                %    break;
                %end
            end
        end
    end
    total_number_of_pulses
    pulses_with_valid_freq
    extracted_average_peak
    size(extracted_average_peak)
    %n = 8; % average every n values
    %b = arrayfun(@(i) mean(extracted_peak(i:i+n-1)),1:n:length(extracted_peak)-n+1)'; % the averaged vector
    %b
            %histogram, mean and SD of peak frequencies
    nbins = 512*1024;
    new_peak=(extracted_average_peak./1000);
    %hold on
    figure
    histogram(new_peak,nbins);
    xlabel('freq(KHz)');
    ylabel('No. of instances');
    title('69KHz fixed freq. tag');
    [Freqs,occurances,ic]=unique(sort(extracted_average_peak));
    Freqs;
    occurances;
    mean_peak_freq=mean(extracted_average_peak)
    sd_freq=std(extracted_average_peak)
end


function [y_fft,y_peak] = fft_pulse(single_pulse,samples,freq_set)
    y_fft=10*log(fftshift(abs(fft(single_pulse))));
    dc_offset=500;
    [M_mag,P_mag]=max(y_fft((samples/2)+dc_offset:samples));
    y_peak=freq_set(samples/2+P_mag+dc_offset);
end
% 
% close all;
% [data,Fs]=audioread(filename);
% time_experiment=180;
% max_index=time_experiment*Fs;
% pulse_length=(10/1000)*Fs;
% threshold=-1/1000;
% threshold_back_samples=0;
% 
% mean(data)
% samples=Fs;
% f=((-samples/2):1:(samples/2)-1)*(Fs/samples);
% 
% figure
% time=0:1/Fs:time_experiment;
% time=time(1:46080000);
% p1=plot(time,data(1:max_index));
% xlabel('time(s)');
% ylabel('Ampl.');
%     %extract first complete envelope
% index=1;
%     while( index < max_index)
%         index=index+1;
%         if data(index) > (threshold)   
%             new_index=index+3.3*Fs; %add 3.3 secs
%             break;
%         end 
%     end
%     %index/Fs
%     %new_index/Fs
%     index=new_index;
%     while( index < max_index)
%         index=index+1;
%         if data(index) > (threshold)   
%             break;
%         end 
%     end
%     %index/Fs
%     %(index-new_index)/Fs
%     if((index-new_index)/Fs>((800/1000)))
%            data(1:index-((1/1000)*Fs))=zeros(index-((1/1000)*Fs),1); 
%     end
% hold on
% p2=plot(time,data(1:max_index));
%     %filterc
% fc = 69000;
% [b,a] = butter(6,fc/(Fs/2));
% dataFilt = filter(b,a,data);
% p3=plot(time,dataFilt(1:max_index));
% h = [p1(1);p2;p3];
% legend(h,'Raw data','Extracted envelopes','Filtered'); %
%     %do fft and other processing
% set=0;
% no_of_sets=500;
% data_extracted=zeros(samples,1,no_of_sets);
% index=1;
%     while( index < max_index)
%         index=index+1;
%         if dataFilt(index) > (threshold)   
%             lower_limit=index-threshold_back_samples;
%             upper_limit=lower_limit+pulse_length-1;
%             set=set+1;
%             if upper_limit > max_index
%                 break;
%             end
%             data_extracted(1:pulse_length,1,set)=dataFilt(lower_limit:upper_limit);
%             index=index+3*pulse_length;
%             %diff=lower_limit-upper_limit;
%             if set >no_of_sets
%                 break;
%             end
%             time_for_one_pulse=0:1/Fs:1;
%             time_for_one_pulse=(lower_limit*(1/Fs))+time_for_one_pulse(1:Fs);
%             hold on
%             %figure
%             plot(time_for_one_pulse, data_extracted(:,1,set),'g');
%         end
%     end
%     set
% fft_set=zeros(samples,1,no_of_sets);
% phase_set=zeros(samples,1,no_of_sets);  
% temp_array=zeros(samples,1);
% %figure
%     for index= 1 : no_of_sets
%         fft_set(:,1,index)=10*log(fftshift(abs(fft(data_extracted(:,1,index)))));
%        % phase_set(:,1,index)=unwrap((angle(fftshift((fft(data_extracted(:,1,index))))))/pi);
%        % threshold = max(abs(phase_set(:,1,index)));
%        % threshold=threshold/10000;
%        % temp_array(abs(phase_set(:,1,index))<threshold) = 0;
%        % phase_set(:,1,index)=temp_array;
%        % plot(f,phase_set(:,1,index));
%        % hold on
%        % figure
%     end
% peak_f=zeros(no_of_sets,1);
%     for index=1:no_of_sets
%         [M_mag,P_mag]=max(fft_set((samples/2)+100:samples,1,index));
%         peak_f(index,1)=f(samples/2+P_mag+100);
%     end
%     peak_f;
%     %[Freqs,occurances,ic]=unique(sort(peak_f));
%    % Freqs;
%    % extracted_peak=zeros(no_of_sets,1);
%     var=0;
%     for index=1:no_of_sets
%         [M_mag,P_mag]=max(fft_set((samples/2):samples,1,index));
%         if(peak_f(index,1)>65880)%68930
%             var=var+1;
%             extracted_peak_old(var,1)=peak_f(index,1);
%         end
%     end
%     var
%     extracted_peak_old;
%     n = 8; % average every n values
%     b = arrayfun(@(i) mean(extracted_peak_old(i:i+n-1)),1:n:length(extracted_peak_old)-n+1)'; % the averaged vector
%     extracted_peak=b
%     size(extracted_peak)
%     nbins = 512*1024;
%     new_peak=(extracted_peak./1000);
%     %hold on
%     figure
%     histogram(new_peak,nbins);
%     xlabel('freq(KHz)');
%     ylabel('No. of instances');
%     title('69KHz fixed freq. tag');
%     [Freqs,occurances,ic]=unique(sort(extracted_peak));
%     Freqs;
%     occurances;
%     mean_peak_freq=mean(extracted_peak)
%     sd_freq=std(extracted_peak)
% end