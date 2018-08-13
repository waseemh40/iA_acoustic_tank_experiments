function z=fft_comparison_filter_v3(filename,experiment_time)
%version 3
%1. Uses 3x pulse length samples inside loop
%Reason is: In moving tag, noise threshold is high and with 1-pulse lenght,
%actual pulse is missed a lot.
%2. Also changed echo rejection time from 300 to 800 ms since this is the case
%for function generator.
close all;
    %raw data extraction
[data,Fs]=audioread(filename);
    %gp variables and parameters
time_experiment=experiment_time*60;
max_index=time_experiment*Fs;
pulse_length=(150/1000)*Fs; %(need more time since noise is above threshold in case of moving tag. Static is OK with exact time/pulse lentgh since noise is lower than thold. )

threshold=.75/1000;     %ideally it should be 1mV but to exclude reflections 2.5mV makes more sense
threshold_back_samples=0;

samples=Fs;
freq_set=((-samples/2):1:(samples/2)-1)*(Fs/samples);

total_number_of_pulses=0;
pulses_with_valid_freq=0;
last_index=(-300*Fs/1000);  %in order to avoid missing first pulse
data_extracted=zeros(samples,1);
single_pulse=zeros(samples,1);
%fft_set=zeros(samples,1,no_of_sets);
%fft_temp=zeros(samples,1);
    %raw data plot
mean_value=mean(data);
data=data-mean_value;
figure
time=0:1/Fs:time_experiment;
time=time(1:end-1);
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
        if dataFilt(index) > (threshold)   
            lower_limit=index-threshold_back_samples;
            if(lower_limit<0)
               index=index+3*pulse_length;
            else
                total_number_of_pulses=total_number_of_pulses+1;
                upper_limit=lower_limit+(3*pulse_length)-1;
                if upper_limit > max_index
                    break;
                end
                data_extracted(1:3*pulse_length,1)=dataFilt(lower_limit:upper_limit);
                single_pulse(:,1)=data_extracted(:,1);
                [fft_temp,temp_f_peak]=fft_pulse(single_pulse,samples,freq_set);
                if temp_f_peak > 66500
                    if((index-last_index)/Fs>((800/1000)))
                    %%
                        pulses_with_valid_freq=pulses_with_valid_freq+1;
                        %fft_set=fft_temp(samples,1,pulses_with_valid_freq);
                        extracted_peak(pulses_with_valid_freq,1)=temp_f_peak;
                            %plot only valid pulses, this will help in next
                            %version of program/algo
                        time_for_one_pulse=0:1/Fs:1;
                        time_for_one_pulse=(lower_limit*(1/Fs))+time_for_one_pulse(1:256000);
                        hold on
                        plot(time_for_one_pulse, data_extracted(:,1),'g');
                        last_index=index;       %used for rejections of echos....
                        index=index+3*pulse_length;
                    else
                        index=index+3*pulse_length; %in case of echo
                    end
                else
                    index=index+(1/1000)*Fs;        %1msec in case of invalid pulse
                end
            end
        end
    end
    total_number_of_pulses
    pulses_with_valid_freq
    extracted_peak;
        %histogram, mean and SD of peak frequencies
    nbins = 512*1024;
    new_peak=(extracted_peak./1000);
    %hold on
    figure
    histogram(new_peak,nbins);
    xlabel('freq(KHz)');
    ylabel('No. of instances');
    title('69KHz fixed freq. tag');
    [Freqs,occurances,ic]=unique(sort(extracted_peak));
    Freqs
    occurances
    mean_peak_freq=mean(extracted_peak)
    sd_freq=std(extracted_peak)
end


function [y_fft,y_peak] = fft_pulse(single_pulse,samples,freq_set)
    y_fft=10*log(fftshift(abs(fft(single_pulse))));
    dc_offset=500;
    [M_mag,P_mag]=max(y_fft((samples/2)+dc_offset:samples));
    y_peak=freq_set(samples/2+P_mag+dc_offset);
end