clear all
close all
clc

%% Loader data 
load('N59_rest.mat');           % Loader EKG- og SKG-data fra forsøgspersonen

tid_original = b1(:,1);         % Definer første kolonne i data, som er tiden
skgAC_original = b1(:,2);       % Definer anden kolonne i data, som er SKG (AC) 
ekg2_original = b1(:,5);        % Definer femte kolonne i data, som er EKG-II
Fs_original = 5000;             % Samplerate som data er optaget med
down_sample_coeffecient = 5;    % Denne koefficent bruges til at downsample EKG- og SKG-data til 1000 Hz.
                              

%% Tilpasser data %% 
data_ekg = b1((2*10^5):(length(b1)-(2*10^5)),5);        % Skærer 200000 samples fra hver ende i EKG-data
data_skg = b1((2*10^5):(length(b1)-(2*10^5)),2);        % Skærer 200000 samples fra hver ende i SKG-data
data_tid = b1((2*10^5):(length(b1)-(2*10^5)),1);        % Skærer 200000 samples fra hver ende i tid-data                                                       % s? arrayl?ngden er den samme
data_skg = (-0.0176*1.1).*data_skg;                     % Ganger en faktor på SKG-signalet, så størrelsesorden 
                                                        % stemmeroverens med det reele/oprendelige SKG-signal
data_skg = data_skg(1:down_sample_coeffecient:end,1);   % Downsampler SKG-signalet til 1000 Hz, ved at beholde hvert 5. datapunkt
data_ekg = data_ekg(1:down_sample_coeffecient:end,1);   % Downsampler EKG-signalet til 1000 Hz, ved at beholde hvert 5. datapunkt
data_tid = data_tid(1:down_sample_coeffecient:end,1);   % Downsampler tid til 1000 Hz, ved at beholde hvert 5. datapunkt

%% Filtrering af EKG   
data_ekg2 = highpass(data_ekg,0.9,Fs_original/down_sample_coeffecient);         % Bruges til at bestemme filterkoeffecienterne A og B
y = zeros(size(data_ekg));      % Definer en vektor af 0 med samme længe som EKG-dataet(downsamlet).

a = [1 -2.98869028151567 2.97744442748529 -0.988753966159255];                  % Filterkoeffecient (udregnet i linje 26)
b = [0.994361084395028 -2.98308325318508 2.98308325318508 -0.994361084395028];  % Filterkoeffecient (udregnet i linje 26)

for nn = 1:length(data_ekg)         % EKG-signalet filteres ud fra filterkoeficenterne via et IIR-filter
                                    % I IIR-filter beregnes outputtet ved bruge af tidligere værdier for output-signalet samtidigtmed værdier fra input-signalet.
    if nn == 1                      % For den første sample bruges kun det nuværende input (enkeltstående tilfælde)
        y(nn) = b(1)*data_ekg(nn);

    elseif nn < numel(b)            % For samples efter den førte (for resten af dataet efter det enkeltstående tilfælde)
        y(nn) = b(1)*data_ekg(nn);  % Udregnet respons for det aktuelle input
            
        for mm = 2:nn               % Beregner responsen ud fra det tidligere indput og output samt filterkoefficenterne.
            y(nn) = y(nn) + b(mm)*data_ekg(nn-mm+1) - a(mm)*y(nn-mm+1); %  Den generelle differensligning for IIR-filtre
        end
        
    % Den generalle filtering, stort set svarende til den overstående filtering. Filterer hele signalet.
    else
        y(nn) = b(1)*data_ekg(nn);  % Udregnet respons for det aktuelle input
    
        for mm = 2:numel(b)         % Beregner responsen ud fra det tidligere indput og output samt filterkoefficenterne.
            y(nn) = y(nn) + b(mm)*data_ekg(nn-mm+1) - a(mm)*y(nn-mm+1);
        end
    end
end

%% Finder R-takker
intsize = 6000/down_sample_coeffecient;                          % Definerer et størrelse 1200 (Bruges til at dele signalet op i segmenter)
y = intsize;
antal_int = floor(length(data_ekg)/intsize)-1;                   % Afrunder længden af EKG-dataet delt med segmentstørrelsen til heltal

for i = 1:antal_int-1                                            % for-loop til identifikation af R-takker og tilhørende x-værdi
 [R_y(i), R_x(i)] = max(data_ekg((1+y):(intsize+y),1));          % Finder x- og y-værdien for R-takken i signalsegmentet på 1200
 R_x(i)=R_x(i)+y;
 y = y+intsize;
end 

%% Finder t-bølger
for i = 1:antal_int-1           % Finder T-bølger ved at finde maks i intervallet fra R-takkens x-værdi + 33 sampels frem til R-takkens x-værdi og 400 sampels frem  
   [t_y(i),t_x(i)] = max(data_ekg(R_x(i)+floor(165/down_sample_coeffecient):R_x(i)+2000/down_sample_coeffecient,1));
   t_x(i) = t_x(i)+R_x(i)+floor(165/down_sample_coeffecient);       % Bestemmer x-værdien til T-bølgens maks
end  

%% For-loop til at aligne signalerne ved brug af xcorr. 
% Dette loop sammenligner kun de to første AC-amplituder
% xcorr sammenligner ved at forskyde signalerne et sample af gangen og finde korrelationsværdien hver gang
aligned_intervals = zeros(length(data_skg(t_x(i):t_x(i)+1200/down_sample_coeffecient)),antal_int-3);

%% For-loop til at aligne signalerne ved brug af xcorr
% Dette for-loop gælder for de resterende intervaller fra 3. og frem
round = 0;             
count = 0;

% Looper et signalsegment(på 10 sampels) gennem alle segmenterne for at bestemme korrelationsværdien 
for k = 1:10:antal_int
    if round == 0    
       
    int1 = data_skg(t_x(k):t_x(k)+1200/down_sample_coeffecient);    % Definerer det første signalsegment af SKG-signalet
    aligned_intervals(:,1) = int1;

    for i = 2:antal_int-1
    int2 = data_skg(t_x(i):t_x(i)+1200/down_sample_coeffecient);    % Definerer det andet signalsegment af SKG-signalet

        [C1,lag1] = xcorr(int1,int2);                               % Returnerer de forsinkelser som bruges til at beregne korrelationen  
        [max_corr_y,max_corr_x] = max(C1);                          % Returnerer index til maks værdien af C1 
        t1 = lag1(max_corr_x);

        if t1 > 0 % Kører hvis intervallet skal forskydes mod højre for at passe med de foregående
            B=zeros(size(int2));
            n= abs(t1); 
            B(n+1:end)=int2(1:end-n);

            aligned_intervals(:,i) = B;

        else % Kører hvis intervallet skal forskydes mod venstre for at passe med de foregående
            B=zeros(size(int1));
            n= abs(t1); 
            B(n+1:end)=int1(1:end-n);  

            aligned_intervals(:,i) = int2;
        end 
    end    

    segment_aligned = corr(aligned_intervals);
    mean_segment_aligned = mean(segment_aligned);
    segmenter = length(mean_segment_aligned);

        for o = segmenter:-1:2                      % Loopet frafilterer signaler, som har en korrelationsfaktor op under 0,7.
           if mean_segment_aligned(1,o) < 0.7;
               aligned_intervals(:,o) = [];
               o = o-1;
           end  
        end

sum = 0;

holder = size(aligned_intervals);
length_aligned_intervals = holder(1,2);

    for i = 1:length_aligned_intervals
    sum = sum + aligned_intervals(:,i);         % Loopet summer de aligned signaler op     
    end

end

    if length_aligned_intervals >= 20           % Når der er taget summen af 20 aligned signaler stopper funktionen 
    round = 1   
    count = count+1
    end
end

    if length_aligned_intervals < 20            % Hvis algoritmen ikke finder 20 aligned signaler, bliver der displayet en fejl
        disp('Error: Not enough ac_ampltiude intervals were detected. Please restart recording') 
        return
    end



sum_mean = sum/(length_aligned_intervals);          % Bestemmer den gennemsnitlige aligned signalsegment 

subplot(3,2,4)                                      % Den gennemsnitlige aligned signalsegment plottes 
plot(sum_mean);                                     

s_range = 100/down_sample_coeffecient               % Definerer et interval på 20 samples

% ac amplituden findes
[max_ac_y,max_ac_x] = max(sum_mean(1:800/down_sample_coeffecient));     % AC-max findes ved at søge efter den maksimale værdi i det mean aligned  
[min_ac_y,min_ac_x] = min(sum_mean((max_ac_x-s_range):max_ac_x,1));     % AC-min findes ved at lede efter den mindste værdi fra AC-maks og 20 samples tilbage 

min_ac_x = min_ac_x+(max_ac_x-s_range);        % AC-min x-værdi bestemmes
ac_amp = abs((min_ac_y-max_ac_y))              % AC-amplituden bestemmes 

%% Overstående algoritme plottes i subplot
hold on
scatter(min_ac_x,min_ac_y);
hold on
scatter(max_ac_x,max_ac_y);
title('Summerede SKG-intervaller(AC-amplitude)')
xlabel('Samples')
ylabel('Amplitude(g)')


subplot(3,2,5)
for i = 10:antal_int-2  
plot(data_skg((t_x(i):t_x(i)+1200/down_sample_coeffecient)));
hold on
end
title('10 ikke-alignede signaler')
xlabel('Samples')
ylabel('Amplitude(g)')

nuller = zeros(1,length(t_x));

subplot(3,2,2)
plot(data_tid,data_skg)
hold on
scatter(data_tid(R_x),nuller) 
hold on
scatter(data_tid(t_x),nuller)
title('SKG med tidspunkt for R-takker og T-bølger')
xlabel('Tid(s)')
ylabel('Amplitude(g)')


subplot(3,2,6)
plot(lag1,C1)
hold on
scatter(t1,max_corr_y)
title('Korrelationsfaktorer')

subplot(3,2,3)
plot(aligned_intervals);
hold on
title('Alignede SKG-intervaller')
xlabel('Samples')
ylabel('Amplitude(g)')

%% Plotter ekg-data med cirkler på R-takker og T-bølger
subplot(3,2,1)
plot(data_tid,data_ekg)            
hold on
scatter(data_tid(R_x),R_y)         %Vi plotter cirkler p? r-takkerne som verficering 
hold on
scatter(data_tid(t_x),t_y);

title('EKG-II')             %titel p? plot
xlabel('Tid(s)')                    %x-akse-titel
ylabel('EKG-II(?)')                 %y-akse-titel     