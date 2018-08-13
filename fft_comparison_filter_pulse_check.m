function z=fft_comparison_filter_pulse_check(filename)

close all;
    %raw data extraction
[data,Fs]=audioread(filename);
    %gp variables and parameters
time_experiment=180;
max_index=time_experiment*Fs;
pulse_length=(10/1000)*Fs;

threshold=-1/1000;
threshold_back_samples=0;

samples=Fs;
f=((-samples/2):1:(samples/2)-1)*(Fs/samples);

total_number_of_pulses=0;
no_of_sets=500;
data_extracted=zeros(samples,1,no_of_sets);
    %raw data plot
figure
time=0:1/Fs:time_experiment;
time=time(1:46080000);
p1=plot(time,data(1:max_index));
xlabel('time(s)');
ylabel('Ampl.');
    %filter and filtered data plot
fc = 69000;
[b,a] = butter(6,fc/(Fs/2));
dataFilt = filter(b,a,data);
hold on
p2=plot(time,dataFilt(1:max_index));
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
                if(index-1+3.5*Fs>max_index)
                    break;
                end
                single_envelope(:,1)=dataFilt(index:index-1+3.5*Fs,1);
                valid_pulse_check=check_pulse(single_envelope,threshold);
                if(valid_pulse_check)
                    total_number_of_pulses=total_number_of_pulses+1;
                    upper_limit=lower_limit+pulse_length-1;
                    if upper_limit > max_index
                        break;
                    end
                    data_extracted(1:pulse_length,1,total_number_of_pulses)=dataFilt(lower_limit:upper_limit);
                    index=index+3*pulse_length;
                    if total_number_of_pulses >no_of_sets
                        break;
                    end
                    time_for_one_pulse=0:1/Fs:1;
                    time_for_one_pulse=(lower_limit*(1/Fs))+time_for_one_pulse(1:256000);
                    hold on
                    plot(time_for_one_pulse, data_extracted(:,1,total_number_of_pulses),'g');
                else
                    index=index+3*pulse_length;
                end
            end
        end
    end
    total_number_of_pulses
        %fft
fft_set=zeros(samples,1,no_of_sets);
%phase_set=zeros(samples,1,no_of_sets);  
%temp_array=zeros(samples,1);
%figure
    for index= 1 : no_of_sets
        fft_set(:,1,index)=10*log(fftshift(abs(fft(data_extracted(:,1,index)))));
       % phase_set(:,1,index)=unwrap((angle(fftshift((fft(data_extracted(:,1,index))))))/pi);
       % threshold = max(abs(phase_set(:,1,index)));
       % threshold=threshold/10000;
       % temp_array(abs(phase_set(:,1,index))<threshold) = 0;
       % phase_set(:,1,index)=temp_array;
       % plot(f,phase_set(:,1,index));
       % hold on
       % figure
    end
peak_f=zeros(no_of_sets,1);
    %peak freq extraction
dc_offset=500;
    for index=1:no_of_sets
        [M_mag,P_mag]=max(fft_set((samples/2)+dc_offset:samples,1,index));
        peak_f(index,1)=f(samples/2+P_mag+dc_offset);
    end
    peak_f;
    [Freqs,occurances,ic]=unique(sort(peak_f));
    Freqs;
    %extracting peaks greater than cut off
cut_off_freq=65880;
    valid_peaks=0;
    for index=1:no_of_sets
        [M_mag,P_mag]=max(fft_set((samples/2):samples,1,index));
        if(peak_f(index,1)>cut_off_freq)
            valid_peaks=valid_peaks+1;
            extracted_peak(valid_peaks,1)=peak_f(index,1);
        end
    end
    valid_peaks
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
    Freqs;
    occurances;
    mean_peak_freq=mean(extracted_peak)
    sd_freq=std(extracted_peak)
end

function y = check_pulse(single_pulse,threshold)
    index=1;
    Fs=256000;
%    time=0:1/Fs:3.5;
%    figure
%    plot(time(1:896000),single_pulse,'g')

    max_index=(10/1000)*Fs;
    instance_above_th=1;
    
    while( index < max_index)
        index=index+1;
        if single_pulse(index) > (threshold)   
            instance_above_th=instance_above_th+1;
        end 
    end
    instance_above_th;
    if((instance_above_th/Fs)>(3/1000))
        y=1;
    else
        y=0;
    end
end