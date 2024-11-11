clear all;
close all;


%% Hard Coded Parameters

% Pass Band parameters

fs = 2000; % Sampling frequency
param_slide = 500; %Window size, peak detection
fb_low = 5; % band-pass low frequency
fb_high = 1000; % band-pass high frequency
ob = 3; % band-pass order

% Notch parameters

fc_notch = 50; % frequency notch filter
q = 10 ; % notch quality factor


%% Muscles Color Representation / Nerve correspondance

% Radial Nerve : ECR, TRI ; 
% Median Nerve: FCR, APB ; 
% Ulnar-Median Nerve : FDP
% Ulnar Nerve : FDI ;
% Axillary Nerve : DELT ;
% Musculocutaneous Nerve : BI ;

% Ulnar-Median Nerve

bleuStart = [0.7, 0.9, 1];  % Light Blue
bleuEnd = [0, 0, 0.4];      % Dark Blue
nbColorsBlue = 4;         % Number of Colors
gradientBlue = [linspace(bleuStart(1), bleuEnd(1), nbColorsBlue)', linspace(bleuStart(2), bleuEnd(2), nbColorsBlue)', linspace(bleuStart(3), bleuEnd(3), nbColorsBlue)'];

% Nerf Radial

orangeStart = [1, 0.7, 0];  % Light Orange
orangeEnd = [1, 0.5, 0];    % Dark Orange
nbColorsOrange = 2;       % Nomber of colors
gradientOrange = [linspace(orangeStart(1), orangeEnd(1), nbColorsOrange)',linspace(orangeStart(2), orangeEnd(2), nbColorsOrange)',linspace(orangeStart(3), orangeEnd(3), nbColorsOrange)'];

% Nerf Auxiliaire

gradientGreen = [0.2, 0.8, 0.2]; % Green

% Nerf Musculo-cutané

gradientBrown = [0.6, 0.3, 0.1]; % Brown

% Trigger

gradientTrigger = [0.502,0.502,0.502];


colors_muscles = struct('ECR',gradientOrange(1,:),'TRI',gradientOrange(2,:),'FCR',gradientBlue(1,:),'ABP',gradientBlue(2,:),'FDP',gradientBlue(3,:),'FDI',gradientBlue(4,:),'DELT',gradientGreen(1,:),'BI',gradientBrown(1,:),'Trigger',gradientTrigger(1,:));







%%%%%%%%%%%




%% Parameters for the user

% Do you want to save the :
save_stim = false; % Stims Bar Graph ?
save_emg = false; % EMG graphs ?
save_area = false; % Area Bar Graph ?

% Do you want to display the :
rep_stim = false; % Stims Bar Graph ?
rep_area = false; % Area Bar Graph ?
rep_test_sliding = false; % Sliding Windows Values ?

% Do you want to apply a :
notch = false; % Notch Filter
pass_band = false; % Band Pass Filter

% Indicate the path containing the recordings files :
folderPath = 'C:\Users\nabil\Desktop\Zahra\FRM labview dvt\NB_20241029_c';

% Indicate the number of the Reference Condition (Task Without Stimulation)
ref_number = 16; 




%%%%%%%%%%%






% Retrieve the different files corresponding to burst experiments
fileNames = {dir(fullfile(folderPath, '*.tdms')).name}; 
indice_test = cell(1, length(fileNames)); % List of condition numbers
for i = 1:length(fileNames)
    condition = regexp(fileNames{i}, 'C(\d+)', 'tokens');
    if ~isempty(condition)
        indice_test{i} = condition{1}{1};
    end
end

% Re-Arrange Conditions List as to have Reference Condition First
ref_index = find(strcmp(indice_test, num2str(ref_number)));
indice_test = [indice_test(ref_index), indice_test([1:ref_index-1, ref_index+1:end])];

param_f_a = cell(1,length(indice_test)); % Test Parameters for each condition : Condition Number / Stimulation Frequency / Stimulation Amplitude


is_ref = false;
for test_num_indice = 1:length(indice_test)
    test_num = indice_test{test_num_indice};

    %% Recollection of TDMS data / Shortcut save as .mat
    
    shortcut_file = ['test_',test_num,'_continuous.mat'];
    if isfile(shortcut_file)
       data_s = load(shortcut_file); % load .mat file if exists
       tdms = data_s.tdms;

    else
        files = dir(fullfile(folderPath, ['C', test_num, '*.tdms']));
        file = fullfile(folderPath, files(1).name); % load .tdms if first reading/no .mat file
        tdms = TDMS_readTDMSFile(file);
        save(shortcut_file,'tdms'); % create .mat shortcut / faster to read
    end
    
    shortcut_file_text = ['test_',test_num,'_texte_continuous.mat'];
    if isfile(shortcut_file_text)
       data_s_texte = load(shortcut_file_text);
       data_texte = data_s_texte.data_texte;
    else
        file_text = [folderPath,'\C',test_num,'.txt'];
        data_texte = readtable(file_text, 'Delimiter', '\t');
        save(shortcut_file_text, 'data_texte');
    end
    
    ReferencePath = [folderPath,'\reference.csv'];
    if isfile(ReferencePath)
        reference = readcell(ReferencePath);
    end
    
    %% Retrieve Stimulation Parameters
    
    freq_stimulation = (data_texte{1,1});
    amplitude_stimulation = (data_texte{1,3});
    stim_alone_duration = (data_texte{1,2});
    if test_num_indice ~= 1
        param_f_a{test_num_indice} = [indice_test{test_num_indice},' - ',num2str(freq_stimulation),'Hz - ',num2str(amplitude_stimulation),'mA'];
    end


    %% Classification Signal Data

    values = tdms.data;
    values_final = values(2:length(values));
    signals_names = tdms.chanNames{1,1}(2:end);
    muscles_names = tdms.chanNames{1,1}(2:end-1);
    n_signals = length(signals_names);
    n_muscles = length(muscles_names);
    
    if test_num_indice == 1  % If Reference Condition ie First Condition in Re-Arranged List (Task without Stimulation)
        donnees_muscles = cell(n_muscles, 4);
        is_ref = true;
        param_f_a{test_num_indice} = 'Without Stim';
    else
        is_ref = false;
    end


    %% Data Processing : Stim and Area (normalized by the duration of mouvement peaks) for each Muscle and Condition
    
    fig1 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

    EMG_signal = cell(1, n_signals);
    EMG_enveloppe = cell(1, n_signals);
    stims = cell(1,n_muscles);
    stims_norm = cell(1,n_muscles);
    areas = cell(1,n_muscles);
    areas_norm = cell(1,n_muscles);
    peaks = cell(2,n_muscles);
    mean_slidings = cell(1,n_muscles);
    
    for j = 1:n_signals
        j = n_signals - j +1 ;
    
        signal_muscle = [];
        for k = 1:((length(values_final)/19)-1)
            signal_muscle = [signal_muscle,values_final{(k-1)*19+j+2}];
        end
        signal_muscle_filtered = filter_function ( signal_muscle , fb_low, fb_high, fs, ob, fc_notch, q, pass_band, notch);
        EMG_signal{j} = signal_muscle_filtered;
    
        %% A. Detection of Stimulation based on Trigger Signal
        if j == n_signals
    
            trigger = EMG_signal{n_signals};
            position0 = find(trigger == 0); % Value of 0 implies a stimulation peak
            diff_trigger = diff(position0); % Number of samples between 2 peaks of value 0
    
            stimulation_start = [];
            stimulation_stop = [];
            int_pulse = floor(fs/freq_stimulation);
    
            stimulation_start(end+1) = position0(1);
            for i = 2:(length(position0)-1)
               
                % If the frequency of 0 values for 2 consecutive pulses differs from the stimulation frequency, it implies a change of burst
                if diff_trigger(i) ~= int_pulse && (diff_trigger(i)) ~= int_pulse+1 && diff_trigger(i) ~= 1
                    stimulation_start(end+1) = position0(i+1);
                    stimulation_stop(end+1) = position0(i);
                end
                
            end
            stimulation_stop(end+1) = position0(length(position0));       
        end  
        

        %% B. Data Calculation
        stims_muscle = [];
        stims_muscle_norm = [];
        area_muscle = [];
        area_muscle_norm = [];
        enveloppe = {};
        start_peak_muscle = [];
        stop_peak_muscle = [];
        mean_sliding_norms_muscle = [];

        if j ~= n_signals
            for k = 1:length(stimulation_start)

                % Retrieve task window
                % After stimulation start and stimulation alone duration
                signal_window = signal_muscle_filtered(stimulation_start(k)+(stim_alone_duration)*fs : stimulation_stop(k));
                length_window = stimulation_stop(k)-(stimulation_start(k)+(stim_alone_duration)*fs);
                
                enveloppe_window = envelope(signal_window,200,'peak'); % Detect the enveloppe of each movement, filtering stimulation frequency
                enveloppe{end+1} = enveloppe_window;
                
                max_muscle = max(signal_window);
                min_muscle = min(signal_window);
                stim_muscle = max_muscle-min_muscle;
                stims_muscle(end+1) = stim_muscle;
                
                if ~is_ref
                    stim_muscle_norm = stim_muscle/reference{1,j+1};
                    stims_muscle_norm(end+1) = stim_muscle_norm;
                    area_norm = area/reference{2,j+1};
                    area_muscle_norm(end+1) = area_norm;
                end
    
                %% C. Movement Peaks Detection using a Sliding Window
                mean_sliding_norms_task = [];
                
                    
                    already_peak = 0;

                    % Construction of the sliding window
                    nb_sliding_windows = floor(length_window/param_slide);
                    for i = 1:nb_sliding_windows
                        if i*param_slide > length_window
                            sliding_window = enveloppe_window((i-1)*param_slide : length_window);
                           end_slide = length_window ;
                        else
                           sliding_window = enveloppe_window((i-1)*param_slide+1 : i*param_slide);
                           end_slide = i*param_slide;
                        end
              
                        mean_sliding_norm = mean(sliding_window)/max(enveloppe_window);
                        mean_sliding_norms_task(end+1) = mean_sliding_norm;
                        
                       
                    end

                    % To facilitate signal comparisons, all signals
                    % returned to x axis
                    mean_sliding_norms_muscle = [mean_sliding_norms_muscle (mean_sliding_norms_task)-min(mean_sliding_norms_task)];

                    % Detection of peaks by percentage discrimination
                    for i = 1:nb_sliding_windows
                        if mean_sliding_norms_muscle((k-1)*nb_sliding_windows+i) > 0.5 && already_peak == 0
                           already_peak = 1;
                           start_peak_muscle(end+1) = (i-1)*param_slide + stimulation_start(k) + (stim_alone_duration)*fs;
                        end
                        if mean_sliding_norms_muscle((k-1)*nb_sliding_windows+i) < 0.5 && already_peak == 1
                           already_peak = 0;
                           stop_peak_muscle(end+1) = (i-1)*param_slide + stimulation_start(k) + (stim_alone_duration)*fs;
                        end
                    
                end

                area = trapz(abs(enveloppe_window));
                area_muscle(end+1) = area;

            end
            
            %% D. Calculation of movement peaks duration
            if length(start_peak_muscle) == length(stop_peak_muscle)  % Movements peaks well detected
                peaks_duration = sum((stop_peak_muscle - start_peak_muscle)/fs);
            else
                peaks_duration = 3*3*length(stimulation_start); % Expected movement durations if default in detection
            end

            %% E. Classification of data
            stims{j} = mean(stims_muscle); 
            areas{j} = mean(area_muscle)/peaks_duration;
            if ~is_ref
                stims_norm{j} = mean(stims_muscle_norm);
                areas_norm{j} = mean(area_muscle_norm);
            end
            EMG_enveloppe{j} = enveloppe;
            peaks{1,j}= start_peak_muscle;
            peaks{2,j}= stop_peak_muscle;
            mean_slidings{j} = mean_sliding_norms_muscle; 

        end
    
        %% Separate Representation between Left and Right Muscles
        a=0;
        if j < ((n_muscles/2))+1
            a=2*j-1;
        else 
            a=2*(j-((n_muscles/2)));
        end


        select = [9,10,11,12,13]; % Muscles targeted by the studied movement
        %for count = 1:length(select) % Uncomment view only concerned muscles
            %if j == select(count)
                subaxis((n_muscles)/2+1,2,a,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)
                %subaxis(length(select),1,count,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)

                    key_color = strtok(signals_names{j}, '_');
                    plot(signal_muscle_filtered,'Color',colors_muscles.(key_color))
                    hold on;
                    if j ~= n_signals
                        for k = 1:length(stimulation_start)
                            plot(stimulation_start(k)+stim_alone_duration*fs : stimulation_start(k)+stim_alone_duration*fs+length(enveloppe{k})-1,enveloppe{k}, 'Color', 'Red')
                        end
                    end
            %end
        %end
        title(signals_names{j},'Color', colors_muscles.(key_color))
        
        xline(stimulation_start,'r--')
        xline(stimulation_stop, 'g--')
        xline(stimulation_start + (stim_alone_duration*fs),'m--')

        for count = 1:length(select)
            if j == select(count)
                xline(start_peak_muscle,'r--')
                xline(stop_peak_muscle,'g--')
            end
        end
        

        %% Re-classification of data by Muscle for Bars Representation
        if j ~= n_signals
            donnees_muscles{j,1}(test_num_indice) = mean(stims_muscle);
            donnees_muscles{j,2}(test_num_indice) = mean(stims_muscle_norm);
            donnees_muscles{j,3}(test_num_indice) = mean(area_muscle);
            donnees_muscles{j,4}(test_num_indice) = mean(area_muscle_norm);
        end
    end


    % Save EMG graphs
    fichier_save = [folderPath,'\EMG'];
    if save_emg
            saveas(fig1, fullfile(fichier_save, ['EMG_',test_num,'.png']));
    end


    %% To analyze in detail the values retrieved by sliding windows
    if rep_test_sliding
        figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
        for j = 9:13
            plot(mean_slidings{1,j});
            hold on;
        end
    end
    
    
    %% MVC / Reference Values Classification and Saving in .csv
    
    % reference.csv contains Stims & Areas Values
    if is_ref
        mean_stim_ref = [];
        mean_area_ref = [];
        for i = 1:n_muscles
            mean_stim_ref(end+1) = mean(stims{i});
            mean_area_ref(end+1) = mean(areas{i});
        end
    
        reference = cell(3,2);
        reference{1,1} = 'Stim Reference';
        reference{2,1} = 'Area Reference';
        reference{1,2} = mean_stim_ref;
        reference{2,2} = mean_area_ref;
        
        writecell(reference,[folderPath,'\reference.csv'])
    end
end


%%%%%%%%%%%


%% Bars Graphs Representation of Stims & Areas

%Stims

if rep_stim
    fichier_save = [folderPath,'\Stim'];
    for k = 1:n_muscles % One figure per muscle

        fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
        subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)

        X = categorical(param_f_a);
        X = reordercats(X,param_f_a); % Gives the stimulation parameters for each condition in legend

        b = bar(X,donnees_muscles{k,1}); % Bar Graph Function

        % Visual Parameters
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        labels = string(b.YData);
        b.FaceColor = 'flat';
        b.CData(1,:) = [.5 0 .5];
        text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
        xlabel('Condition number');
        ylabel(['Stim',muscles_names(k)]);
    
        if save_stim
            saveas(fig, fullfile(fichier_save, ['Stim_',num2str(k),'.png']));
            saveas(fig, fullfile(fichier_save,['Stim_',muscles_names(k),'.pdf']));
            saveas(fig, fullfile(fichier_save, ['Stim_',muscles_names(k),'.svg']));
            savefig(fig , fullfile(fichier_save, ['Stim_',muscles_names(k),'.fig']));
        end
    end
end


%%%%%


%Area
if rep_area
    fichier_save = [folderPath,'\Area'];
    for k = 1:n_muscles % One figure per muscle

        fig = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
        subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 

        X = categorical(param_f_a); % Gives the stimulation parameters for each condition in legend
        X = reordercats(X,param_f_a);

        b = bar(X,donnees_muscles{k,3}); % Bar Graph Function

        % Visual Parameters
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        labels = string(b.YData);
        b.FaceColor = 'flat';
        b.CData(1,:) = [.5 0 .5];
        text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
        xlabel('Condition number');
        ylabel(['Area',muscles_names(k)]);
    
        if save_area
            saveas(fig, fullfile(fichier_save, ['Area_',num2str(k),'.png']));
            saveas(fig, fullfile(fichier_save,['Area_',muscles_names(k),'.pdf']));
            saveas(fig, fullfile(fichier_save, ['Area_',muscles_names(k),'.svg']));
            savefig(fig , fullfile(fichier_save, ['Area_',muscles_names(k),'.fig']));
        end
    end
end


%%%%%%%%%%%


%% Signal Filter Function

function [V_filtered] = filter_function  ( V_raw , fb_low, fb_high, fs, ob, fc_notch, q, pass_band, notch)


        % Band-pass filter
        if pass_band == true 
            Wn = [fb_low fb_high] / (fs/2);
            [b_butter, a_butter] = butter(ob, Wn, 'bandpass');
            V_bp = filter(b_butter, a_butter, V_raw);
        else
            V_bp = V_raw;
        end

        % Notch filter
        if notch == true
            bw = fc_notch/q;
            [b_notch, a_notch] = iirnotch(fc_notch/(fs/2), bw/(fs/2));
            V_filtered = filter(b_notch, a_notch, V_bp);
        else
            V_filtered = V_bp;
        end
end


