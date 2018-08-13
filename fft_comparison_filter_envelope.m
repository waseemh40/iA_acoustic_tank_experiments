function z=fft_comparison_filter_envelope(filename)

close all;
[data,Fs]=audioread(filename);
time_experiment=180;
max_index=time_experiment*Fs;
pulse_length=(10/1000)*Fs;
envelope_length=3.5*Fs;
inter_envelope_time=7*Fs;
half_inter_envelope_time=3.5*Fs;
threshold=0.00;
threshold_back_samples=0;

samples=envelope_length;
f=((-samples/2):1:(samples/2)-1)*(Fs/samples);

set=0;
no_of_sets=50;
figure
time=0:1/Fs:time_experiment;
time=time(1:46080000);
p1=plot(time,data(1:max_index));
xlabel('time(s)');
ylabel('Ampl.');
%figure
data_extracted=zeros(samples,1,no_of_sets);
    %filterc
fc = 69000;

[b,a] = butter(6,fc/(Fs/2));
%freqz(b,a);
dataFilt = filter(b,a,data);
%figure
hold on
p2=plot(time,dataFilt(1:max_index));
h = [p1(1);p2];
legend(h,'Raw data','Filtered'); %
index=1;
    while( index < max_index)
        index=index+1;
        if dataFilt(index) > (threshold)   
            lower_limit=index-threshold_back_samples;
            if(lower_limit<0)
               index=index+half_inter_envelope_time;
            else
                set=set+1;
                upper_limit=lower_limit+envelope_length-1;
                if upper_limit > max_index
                    break;
                end
                data_extracted(1:envelope_length,1,set)=dataFilt(lower_limit:upper_limit);
                index=index+half_inter_envelope_time;
                %diff=lower_limit-upper_limit;
                if set >no_of_sets
                    break;
                end
                time_for_one_pulse=0:1/Fs:3.5;
                time_for_one_pulse=(lower_limit*(1/Fs))+time_for_one_pulse(1:envelope_length);
                hold on
                %figure
                plot(time_for_one_pulse, data_extracted(:,1,set));
            end
        end
    end
    set
fft_set=zeros(samples,1,no_of_sets);
phase_set=zeros(samples,1,no_of_sets);  
temp_array=zeros(samples,1);
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
    for index=1:no_of_sets
        [M_mag,P_mag]=max(fft_set((samples/2)+100:samples,1,index));
        peak_f(index,1)=f(samples/2+P_mag+100);
    end
    peak_f;
    [Freqs,occurances,ic]=unique(sort(peak_f));
    Freqs
   % extracted_peak=zeros(no_of_sets,1);
    var=0;
    for index=1:no_of_sets
        [M_mag,P_mag]=max(fft_set((samples/2):samples,1,index));
        if(peak_f(index,1)>65880)%68930
            var=var+1;
            extracted_peak(var,1)=peak_f(index,1);
        end
    end
    var
    extracted_peak;
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
    x=[68949 68949];
    y=[-200 -10];
    %line(x,y);
    x=[68959 68959];
    %line(x,y);
end