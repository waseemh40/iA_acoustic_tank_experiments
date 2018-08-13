function z=fft_comparison(filename)

[data,Fs]=audioread(filename);
samples=256000;
f=((-samples/2):1:(samples/2)-1)*(Fs/samples);
max_index=4.6*10000000;
set=0;
no_of_sets=140;
data_extracted=zeros(samples,1,no_of_sets);
index=1;
    while( index < max_index)
        index=index+1;
        if data(index) > (0.05)   
            set=set+1;
            lower_limit=index-100;
            upper_limit=index+10000-101;;%samples/8->10000;
            data_extracted(1:+10000,1,set)=data(lower_limit:upper_limit);
            index=index+samples/8;
            %index
            if set >no_of_sets
                break;
            end
            %plot( data_extracted(:,1,set));
            hold on
        end
    end
fft_set=zeros(samples,1,no_of_sets);    
    for index= 1 : no_of_sets
        fft_set(:,1,index)=10*log(fftshift(abs(fft(data_extracted(:,1,index)))/samples));
        %plot(f,fft_set(:,1,index));
        hold on
    end
peak_f=zeros(no_of_sets,1);
    for index=1:no_of_sets
        [M_mag,P_mag]=max(fft_set(samples/2:samples,1,index));
        peak_f(index,1)=f(samples/2+P_mag);
    end
    peak_f
    nbins = 512*1024;
    new_peak=(peak_f./1000);
    histogram(new_peak,nbins);
    xlabel('freq(KHz)');
    ylabel('No. of instances');
    title('69KHz fixed freq. tag');
    [Freqs,occurances,ic]=unique(sort(peak_f));
    Freqs
    occurances
    mean_peak_freq=mean(peak_f)
    sd_freq=std(peak_f)
    x=[68949 68949];
    y=[-200 -10];
    %line(x,y);
    x=[68959 68959];
    %line(x,y);
end