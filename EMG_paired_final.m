clear all;
close all;

%% Hard Coded Parameters :

% Acquisition parameters

fs = 20000; % Sampling frequency
t_acq = 0.2; % Acquisition duration for each pulse
t_acq_d = [-0.05 0.15]; % Time frame aligned with first pulse
d_int_p = 0.1; % Duration between two pulses
d_pulse = 0.0005; % Pulse Width
incr = 5; % Amplitude Increment (mA)

% Pass Band parameters

fb_low = 5; % Band-pass Low Frequency
fb_high = 1000; % Band-pass High Frequency
ob = 3; % Band-Pass Order

% Notch parameters

fc_notch = 50; % Notch Filter Frequency
q = 10 ; % Notch Quality Factor

% Stim Calculation Parameter
t_stim = [0.01 0.04]; % Time Frame for Calculations
discrim_diff = 2; % Discrimination Parameter for Post Activation Depression
discrim_diff_norm = 3; % Normalized Discrimination Parameter for Post Activation Depression

% Threshold Parameters
param_thresh = 1;
param_thresh_norm = 1;

% Normalisation Parameters
t_norm = [-0.05 -0.01]; % Time window for Normalisation

% Representation Parameters
t_rep = [-0.02 0.15] ; % Time Frame Represented
lim_emg = [-0.5 0.5] ; % Amplitude Frame EMG
lim_emg_norm = [-120 120]; % Amplitude Frame EMG Normalized
lim_stim = [0 10]; % Amplitude Frame Stim
lim_recruit = [0 1]; % Amplitude Frame Recruitment Curve
lim_recruit_norm = [0 600]; % Amplitude Frame Normalized Recruitment Curve
lim_diff = [-0.3 0.5]; % Amplitude Frame Difference First and Second Pulse Curve
lim_diff_norm = [-100 150]; % Amplitude Frame Normalized Difference First and Second Pulse Curve
amp_range = [5 100]; % Stimulation Amplitude Range Representation

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






%%%%%%%%%%%%%%%%%%




%% Parameters for the user

% Do you want to apply a :
pass_band = false; % Pass-Band Filter
notch = false; % Notch Filter

% Do you want to normalize by either :
b_rms = false ; % RMS
%OR
b_zscore = true ; % ZScore

% Do you want to display the curves ?
plot_check = true;

% Do you want to save the :
save_emg = false ; % EMG Curve ?
save_emg_norm = false ; % Normalized EMG Curve ?
save_recruit = false ; % Recruitment Curve ?
save_recruit_norm = false ; % Normalized Recruitment Curve ?
save_diff = false ; % Difference First and Second Pulse Curve ?
save_diff_norm = false ; % Normalized Difference First and Second Pulse Curve ?
save_tab = false ; % Threshold Tab ?
save_tab_norm = false ; % Normalized Threshold Tab ?

% Do you want to add a legend for muscles (Recruit / Diff)
legende_m = false;

% Do you want to process and display Post Activation Depression Curves ?
specific = false;
specific_norm = false; % Normalized

% Indicate the path containing the recordings files and the condition number :
folderPath = 'C:\Users\nabil\Desktop\Zahra\FRM labview dvt\NB_20241029_p - Copy\';
test_num = '1';

% Specify in the right order the muscle names
nom_muscles = {'FCR_L','FDP_L','ECR_L','FDI_L','ABP_L','BI_L','TRI_L','DELT_L','FCR_R','FDP_R','ECR_R','FDI_R','ABP_R','BI_R','TRI_R','DELT_R'};

% Specify the list of tested conditions and current condition
condition_list = {'C7/T1','C5/C6','C7/T1 Normalisé','C5/C6 Normalisé'};
num_condition = 2;




%%%%%%%%%%%%%%%%%%







%% Data Recuperation and Creation of Shortcut .mat



shortcut_file = ['test_',test_num,'.mat'];
if isfile(shortcut_file)
   data_s = load(shortcut_file); % load .mat file if exists
   data= data_s.data;
else
    file = [folderPath,test_num,'.txt'];
    opts = detectImportOptions(file, 'Delimiter', '\t'); % load file if first reading / no .mat file
    opts.DataLines = [23 Inf];
    opts.VariableNamesLine = 1;
    data = readtable(file,opts);
    save(shortcut_file,'data'); % create .mat shortcut / faster to read
end


shortcut_file_texte = ['test_',test_num,'_texte.mat'];
if isfile(shortcut_file_texte)
   data_s_texte = load(shortcut_file_texte);
   data_texte = data_s_texte.data_texte;
else
    file_texte = [folderPath,test_num,'txt'];
    data_texte = readtable(file_texte, 'Delimiter', '\t');
    save(shortcut_file_texte, 'data_texte');
end

file_save = [folderPath,test_num];



%% Stimulation Amplitude Parameters

range_impulsion= data_texte{:,2}; % Stimulation Amplitudes
pulse = data_texte{:,3}; % Pulses Index

val0 = range_impulsion(1); % Minimum Amplitude (mA)
amp_max = range_impulsion(length(range_impulsion)); % Maximum Amplitude (mA)
nb_pulse = length(pulse); % Number of Paired-Pulse Stimulation during Acquisition
nb_valeurs = length(unique(range_impulsion)); % Nomber of different Stimulation Amplitude Values



%% Data Classification
n_muscles = length(nom_muscles);
EMG_cell = cell(1,n_muscles);


%%%%%%%%%%%


%% EMG Curve / Overlay for each Stimulation Amplitude

fig1 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

for i = 1:n_muscles

    signal_raw = data{2:end,i+1};
    signal_filtered = filter_function (signal_raw, fb_low, fb_high, fs, ob, fc_notch, q, pass_band, notch);

    % Reshape the signal to separate by paired_pulses
    EMG_cell{i} = reshape(signal_filtered(1:nb_pulse*t_acq*fs),[t_acq*fs,nb_pulse]);

    % Time Vector Construction
    x = linspace (t_acq_d(1)*1000,t_acq_d(2)*1000,t_acq*fs);  

    % Classification Left or Right muscle
    j=0;
    p=0;
    if i < (n_muscles/2)+1
        j=2*i-1;
        p=i;
    else 
        j=2*(i-(n_muscles/2));
        p=i-(n_muscles/2);
    end
    
    % EMG Representation Plot
    subaxis(n_muscles/2,2,j,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)
    colormap_muscle = jet(length(EMG_cell{i}(1,:)));
    
    for k = 1:length(EMG_cell{i}(1,:))
        if plot_check
            plot(x,EMG_cell{i}(:,k),'Color',colormap_muscle(k,:))
            hold on;
        end
    end

    if i == n_muscles/2
        xlabel('Time (ms)');
        ylabel('EMG Response (mV)')
    end

    ylim(lim_emg);
    xlim(t_rep*1000)
    xline(0,'r--'); % Start First Pulse
    xline((d_pulse)*1000,'r--') % End First Pulse
    xline(t_stim(1)*1000,'b--'); % Start Analysis Window First Pulse
    xline(t_stim(2)*1000,'b--'); % End Analysis Window Second Pulse
    xline((d_int_p)*1000,'r--'); % Start Second Pulse
    xline((d_int_p+d_pulse)*1000,'r--'); % End Second Pulse
    xline((t_stim(1)+d_int_p)*1000,'b--'); % Start Analysis Window Second Pulse
    xline((t_stim(2)+d_int_p)*1000,'b--'); % End Analysis Window Second Pulse

    key_color = strtok(nom_muscles{i}, '_');
    title(nom_muscles{i},'Color',colors_muscles.(key_color))
    box off;
end

if save_emg
    saveas(fig1, fullfile(file_save, 'time.png'));
    saveas(fig1, fullfile(file_save, 'time.pdf'));
    saveas(fig1, fullfile(file_save, 'time.svg'));
    savefig(fig1, fullfile(file_save, 'time.fig'));
end


%%%%%%%%%%%


%% Normalized EMG Curve / Overlay for each Stimulation Amplitude

if b_rms || b_zscore % Only if Normalization wanted

    EMG_cell_norm = cell(1,n_muscles);
    
    fig2 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    for i = 1:n_muscles 
    
        EMG_cell_norm{i}=[];
        
        for k = 1:length(EMG_cell{i}(1,:))
    
            signal = EMG_cell{i} (:,k);
            pt_norm = [round((t_norm(1)- t_acq_d(1))*fs) round((t_norm(2)- t_acq_d(1))*fs)]; % Index of baseline window without stim
            baseline = signal(pt_norm(1)+1 : pt_norm(2)+1); % Signal Window without Stimulation
    
            if b_rms % Normalization by RMS
                rms_valeur = rms(baseline) ;
                normalize = signal / rms_valeur;
            end    
    
            if b_zscore % Normalization by ZScore
                normalize = (signal-mean(baseline))/std(baseline);
            end
    
            
            EMG_cell_norm{i} (:,end+1)= normalize;
    
    
            % Classification Left and Right Muscles
            j=0;
            p=0;
            if i < (n_muscles/2)+1
                j=2*i-1;
                p=i;
            else 
                j=2*(i-(n_muscles/2));
                p=i-(n_muscles/2);
            end
            
            % Normalized EMG Representation Plot
            subaxis(n_muscles/2,2,j,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)
            colormap_muscle = jet(length(EMG_cell{i}(1,:)));
             if plot_check
                plot(x,EMG_cell_norm{i}(:,k),'Color',colormap_muscle(k,:))
                hold on;
             end
    
        end
    
        if i == n_muscles/2
                xlabel('Time (ms)');
                ylabel('Normalized EMG Response (.)')
        end
    
        ylim(lim_emg_norm)
        xlim(t_rep*1000)
        xline(0,'r--'); % Start First Pulse
        xline((d_pulse)*1000,'r--') % End First Pulse
        xline(t_stim(1)*1000,'b--'); % Start Analysis Window First Pulse
        xline(t_stim(2)*1000,'b--'); % End Analysis Window Second Pulse
        xline((d_int_p)*1000,'r--'); % Start Second Pulse
        xline((d_int_p+d_pulse)*1000,'r--'); % End Second Pulse
        xline((t_stim(1)+d_int_p)*1000,'b--'); % Start Analysis Window Second Pulse
        xline((t_stim(2)+d_int_p)*1000,'b--'); % End Analysis Window Second Pulse

        key_color = strtok(nom_muscles{i}, '_');
        title(nom_muscles{i},'Color',colors_muscles.(key_color))
        box off;
    end
end

if save_emg_norm
    saveas(fig2, fullfile(file_save, 'time_norm.png'));
    saveas(fig2, fullfile(file_save, 'time_norm.pdf'));
    saveas(fig2, fullfile(file_save, 'time_norm.svg'));
    savefig(fig2 , fullfile(file_save, 'time_norm.fig'));
end


%%%%%%%%%%%


%% Stims Calculations

stims= [];

for i = 1:n_muscles
    stim_muscle = [];
    stim_window = EMG_cell{i}(round((t_stim(1) - t_acq_d(1))*fs):round((t_stim(2) - t_acq_d(1))*fs), 1:nb_pulse); % Analysis window for stims calculation
    min_pulse = min(stim_window);
    max_pulse = max(stim_window);
    stim_pulse = max_pulse - min_pulse;
    stim_muscle = [stim_muscle,stim_pulse] ;
    stims = [stims;stim_muscle] ;
end


%%%%%%%%%%%


%% Normalized Stims Calculations


if b_rms || b_zscore
    
    stims_norm= [];
    
    for i = 1:n_muscles
        stim_muscle_norm = [];
        stim_window_norm = EMG_cell_norm{i}(round((t_stim(1) - t_acq_d(1))*fs):round((t_stim(2) - t_acq_d(1))*fs), 1:nb_pulse); % Analysis window for normalized stims calculation
        min_pulse_norm = min(stim_window_norm);
        max_pulse_norm = max(stim_window_norm);
        stim_pulse_norm = max_pulse_norm - min_pulse_norm;
        stim_muscle_norm = [stim_muscle_norm,stim_pulse_norm] ;
        stims_norm = [stims_norm;stim_muscle_norm] ;
    end

end


%%%%%%%%%%%


%% Recruitment Curve

fig3 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');

thresh = [];
means = [];

for i = 1:n_muscles

   %% Table of Stims Calculations classified by Stimulation Amplitude
   amp_values = unique (range_impulsion); % List of different amplitude values
   stims_arr = [];
   j=1;
   a=0;
   for k = 1:length(stims(1,:))
       if range_impulsion(k) == amp_values(j) % For each Stimulation Amplitude change, the column changes
           a=a+1;
           stims_arr(a,j) = stims(i,k);
       else
                  j=j+1;
                  a=1;
                  stims_arr(a,j) = stims(i,k);
       end
   end

   mean_pulse = mean(stims_arr, 'omitnan'); % Mean of stims over pulses at the same amplitude
   means(:,end+1) = mean_pulse;
    

   %% Calculation of Stimulation Amplitude thresholds of Muscle Activation
   h = 2;
   discrim = [];
   % Discrimation on response amplitude
   while ((mean_pulse(h)-mean_pulse(h-1))/abs(mean_pulse(h-1)))< param_thresh && h<length(mean_pulse) 
       discrim (end+1) = ((mean_pulse(h)-mean_pulse(h-1))/abs(mean_pulse(h-1)));
       h = h+1;
   end
    
   if h<length(mean_pulse) 
       % Discimination on growth of response
       if mean_pulse(h+1) < mean_pulse(h)
          h=h+1;
          while ((mean_pulse(h)-mean_pulse(h-1))/abs(mean_pulse(h-1)))< param_thresh && h<length(mean_pulse)
          discrim (end+1) = ((mean_pulse(h)-mean_pulse(h-1))/abs(mean_pulse(h-1)));
          h = h+1;
          end
       end
   end
   % If the amplitude meets the criteria :
   if h<=length(mean_pulse) && (mean_pulse(length(mean_pulse))-mean_pulse(1))/(mean_pulse(1))>6
       thresh = [thresh amp_values(h)];
   % If no amplitude meets the criteria :
   else
       thresh = [thresh NaN];
   end
   
   key_color = strtok(nom_muscles{i}, '_');
   if i < (n_muscles/2)+1
        if plot_check
           subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07)
           h1 = plot(amp_values,mean_pulse,'--o','DisplayName',nom_muscles{i},'Linewidth',2,'Color',colors_muscles.(key_color));
           % Buttons to control the display of each muscle
           uicontrol('Style', 'checkbox', 'String', ['Afficher ',nom_muscles{i}],'Position', [50, 400+20*i, 150, 20], 'Value', 1, 'Callback', @(src, event) set(h1, 'Visible', boolToOnOff(src.Value)));
           % Display Muscles Names on each curve
           text(amp_values(end) + 0.5, mean_pulse(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
           hold on ;
           % Plot of single points corresponding to stim of each pulse
           scatter(range_impulsion,stims(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
           hold on;
        end
   else
         if plot_check
           subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
           h1 = plot(amp_values,mean_pulse,'-o','DisplayName',nom_muscles{i},'LineWidth',2,'Color',colors_muscles.(key_color));
           % Buttons to control the display of each muscle
           uicontrol('Style', 'checkbox', 'String', ['Afficher ',nom_muscles{i}],'Position', [50, 400+20*i, 150, 20], 'Value', 1, 'Callback', @(src, event) set(h1, 'Visible', boolToOnOff(src.Value)));
           % Display Muscles Names on each curve
           text(amp_values(end) + 0.5, mean_pulse(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
           hold on ;
           % Plot of single points corresponding to stim of each pulse
           scatter(range_impulsion,stims(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
           hold on;
         end
   end

end
    
xlabel('Stimulation intensity (mA)')
xlim(amp_range)
ylim(lim_recruit)
ylabel('Mean EMG response  (mV)')
if legende_m % Controls the display of legend box
    legend
end
title('Recruitment Curve','Color','red')
xticks(val0:incr:amp_max)

if save_recruit
    saveas(fig3, fullfile(file_save, 'recruit.png'));
    saveas(fig3, fullfile(file_save, 'recruit.pdf'));
    saveas(fig3, fullfile(file_save, 'recruit.svg'));
    savefig(fig3 , fullfile(file_save, 'recruit.fig'));
end


%%%%%%%%%%%


%% Threshold Table

rank = round(tiedrank(thresh))-1; % Returns ranks of activation for each muscle

% Representation Parameters Table
table_width = 82 * n_muscles;
table_height = 60;

fig4 = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, table_width, table_height]);

% Table with the activation amplitude and activation rank for each muscle
uitable('Parent', fig4,'Data',[num2cell(thresh); num2cell(rank)],'ColumnName',nom_muscles,'RowName',{'Thresholds','Rank'},'Position', [0 0 table_width table_height]);

if save_tab
    saveas(fig4, fullfile(file_save, 'table.png'));
end


%%%%%%%%%%%


%% Normalized Recruitment Curve

if b_rms || b_zscore

    fig5 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    
    
    thresh_norm = [];
    means_norm = [];
    for i = 1:n_muscles
    
       %% Table of Stims Calculations classified by Stimulation Amplitude
       stims_arr_norm = [];
       j=1;
       a=1;
       for k = 1:length(stims_norm(1,:))
           if range_impulsion(k) == amp_values(j) % For each Stimulation Amplitude change, the column changes
               a=a+1;
               stims_arr_norm(a,j) = stims_norm(i,k);
           else
               j=j+1;
               a=1;
               stims_arr_norm(a,j) = stims_norm(i,k);
           end
       end
    
    
       mean_pulse_norm = mean(stims_arr_norm, 'omitnan'); % Mean of stims over pulses at the same amplitude
       means_norm(:,end+1) = mean_pulse_norm;
       
    
       %% Calculation of Stimulation Amplitude thresholds of Muscle Activation
       h = 2;
       discrim = [];
       % Discrimation on response amplitude
       while ((mean_pulse_norm(h)-mean_pulse_norm(h-1))/abs(mean_pulse_norm(h-1)))< param_thresh_norm && h<length(mean_pulse_norm) 
           discrim (end+1) = ((mean_pulse_norm(h)-mean_pulse_norm(h-1))/abs(mean_pulse_norm(h-1)));
           h = h+1;
       end
        
       if h<length(mean_pulse_norm)
           % Discimination on growth of response
           if mean_pulse_norm(h+1) < mean_pulse_norm(h)
              h=h+1;
              while ((mean_pulse_norm(h)-mean_pulse_norm(h-1))/abs(mean_pulse_norm(h-1)))< param_thresh_norm && h<length(mean_pulse_norm)
              discrim (end+1) = ((mean_pulse_norm(h)-mean_pulse_norm(h-1))/abs(mean_pulse_norm(h-1)));
              h = h+1;
              end
            end
       end
       % If the amplitude meets the criteria :
       if h<=length(mean_pulse_norm) && (mean_pulse_norm(length(mean_pulse_norm))-mean_pulse_norm(1))/(mean_pulse_norm(1))>6 
            thresh_norm = [thresh_norm amp_values(h)];
       % If no amplitude meets the criteria :
       else
           thresh_norm = [thresh_norm NaN];
       end
    
       key_color = strtok(nom_muscles{i}, '_');
           if i < (n_muscles/2)+1
                if plot_check
                   subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
                    h1 = plot(amp_values,mean_pulse_norm,'--o','DisplayName',nom_muscles{i},'Linewidth',2,'Color',colors_muscles.(key_color));
                    % Buttons to control the display of each muscle
                    uicontrol('Style', 'checkbox', 'String', ['Afficher ',nom_muscles{i}],'Position', [50, 400+20*i, 150, 20], 'Value', 1, 'Callback', @(src, event) set(h1, 'Visible', boolToOnOff(src.Value)));
                    % Display Muscles Names on each curve
                    text(amp_values(end) + 0.5, mean_pulse_norm(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
                    hold on ;
                    % Plot of single points corresponding to stim of each pulse
                    scatter(range_impulsion,stims_norm(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
                    hold on;
                end
               
           else
                if plot_check
                    subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
                    h1 = plot(amp_values,mean_pulse_norm,'-o','DisplayName',nom_muscles{i},'LineWidth',2,'Color',colors_muscles.(key_color));
                    % Buttons to control the display of each muscle
                    uicontrol('Style', 'checkbox', 'String', ['Afficher ',nom_muscles{i}],'Position', [50, 400+20*i, 150, 20], 'Value', 1, 'Callback', @(src, event) set(h1, 'Visible', boolToOnOff(src.Value)));
                    % Display Muscles Names on each curve
                    text(amp_values(end) + 0.5, mean_pulse_norm(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
                    hold on ;
                    % Plot of single points corresponding to stim of each pulse
                    scatter(range_impulsion,stims_norm(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
                    hold on;
                end
           end
    
    end
    xlabel('Stimsulation intensity (mA)')
    xlim(amp_range);
    ylim(lim_recruit_norm);
    ylabel('Mean normalized EMG response  (.)')
    if legende_m % Controls the display of legend box
        legend
    end
    title('Recruitment Curve of the normalized signal','Color','red')
    xticks(val0:incr:amp_max)
    
    if save_recruit_norm
        saveas(fig5, fullfile(file_save, 'recruit_norm.png'));
        saveas(fig5, fullfile(file_save, 'recruit_norm.pdf'));
        saveas(fig5, fullfile(file_save, 'recruit_norm.svg'));
        savefig(fig5 , fullfile(file_save, 'recruit_norm.fig'));
    end
end




%%%%%%%%%%%


%% Normalized Threshold Table

if b_rms || b_zscore
    rank_norm = round(tiedrank(thresh_norm))-1;% Returns ranks of activation for each muscle
    
    % Representation Parameters Table
    table_width = 82 * n_muscles;
    table_height = 60;
    
    fig6 = figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, table_width, table_height]);
    
    % Table with the activation amplitude and activation rank for each muscle
    uitable('Parent', fig6 ,'Data',[num2cell(thresh_norm); num2cell(rank_norm)],'ColumnName',nom_muscles,'RowName',{'Thresholds','Rank'},'Position', [0 0 table_width table_height]);
    
    
    if save_tab_norm
        saveas(fig6, fullfile(file_save, 'table_norm.png'));
    end
end

%%%%%%%%%%%


%% Calculations First and Second Activation for each Paired-Pulse

diff = [];
for i = 1:n_muscles
   diff_muscle = [];
   
   for k = 1:nb_pulse
       signal_diff = EMG_cell{i}(:,k);
       pulses_1 = signal_diff(round((-t_acq_d(1)+t_stim(1))*fs):round((-t_acq_d(1)+t_stim(2))*fs)); % Retrieve the signal from the first analysis window
       pulses_2 = signal_diff(round((-t_acq_d(1)+t_stim(1)+d_int_p)*fs):round((-t_acq_d(1)+t_stim(2)+d_int_p)*fs)); % Retrieve the signal from the first analysis window
       stim_1= (max(pulses_1)-min(pulses_1)); % Stim from first analysis window
       stim_2= (max(pulses_2)-min(pulses_2)); % Stim from second analysis window
       diff_p = stim_1 - stim_2; % Calculation of the difference for each paired-pulse (Other posibility: ratio)
       diff_muscle(end+1)=diff_p;
       
   end
   diff(end+1,:)=diff_muscle;
end


%%%%%%%%%%%


%% Post Depression Activation Comparison Curve

fig7 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
diff_arr = {};

for i = 1:n_muscles

   %% Table of Difference/Ratio Calculations classified by Stimulation Amplitude
   j=1;
   a=0;
   for k = 1:length(stims)
        if range_impulsion(k) == amp_values(j) % For each Stimulation Amplitude change, the column changes
           a=a+1;
            diff_muscle_arr(a,j) = diff(i,k);
        else
            j=j+1;
            a=1;
            diff_muscle_arr(a,j) = diff(i,k);
        end
   end
   diff_arr{end+1}=diff_muscle_arr;
   mean_diff = mean(diff_muscle_arr, 'omitnan'); % Mean of differences/ratios over pulses at the same amplitude

   key_color = strtok(nom_muscles{i}, '_');
   if i < (n_muscles/2)+1
        if plot_check
            subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
            plot(amp_values,mean_diff,'--o','DisplayName',nom_muscles{i},'Linewidth',2,'Color',colors_muscles.(key_color))
            % Display Muscles Names on each curve
            text(amp_values(end) + 0.5, mean_diff(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
            hold on ;
            % Plot of single points corresponding to difference/ratio of each pulse
            scatter(range_impulsion,diff(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
            hold on;

        end
   else
        if plot_check
            subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
            plot(amp_values,mean_diff,'-o','DisplayName',nom_muscles{i},'LineWidth',2,'Color',colors_muscles.(key_color))
            % Display Muscles Names on each curve
            text(amp_values(end) + 0.5, mean_diff(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
            hold on ;
            % Plot of single points corresponding to difference/ratio of each pulse
            scatter(range_impulsion,diff(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
            hold on;
        end
   end

end
xlabel('Stimulation intensity (mA)')
xlim(amp_range);
ylim(lim_diff);
ylabel('Mean difference first and second pulse amplitude  (mV)')
if legende_m % Controls the display of legend box
    legend
end
title('Mean difference first and second pulse amplitude over time','Color','red')
xticks(val0:incr:amp_max)



if save_diff
    saveas(fig7, fullfile(file_save, 'diff.png'));
    saveas(fig7, fullfile(file_save, 'diff.pdf'));
    saveas(fig7, fullfile(file_save, 'diff.svg'));
    savefig(fig7 , fullfile(file_save, 'diff.fig'));
end



%%%%%%%%%%%


%% Post Activation Depression Selection Curve

if specific
    
    liste_sel={};
    for i = 1:n_muscles
            
        first_pulse = EMG_cell{i}(round((t_stim(1)-t_acq_d(1))*fs):round((t_stim(2)-t_acq_d(1))*fs), 1:nb_pulse); % Contains the first pulse analysis window of each paired-pulse stimulation
        second_pulse = EMG_cell{i}(round((t_stim(1)+d_int_p-t_acq_d(1))*fs):round((t_stim(2)+d_int_p-t_acq_d(1))*fs), 1:nb_pulse); % Contains the second pulse analysis window of each paired-pulse stimulation
        
        diff_muscle_arr = diff_arr{i};
        %% Selection of paired-pulses with sufficient Post Activation Depression
        liste_sel_muscle = [];
        for k = 1:length(diff_muscle_arr(1,:))
            
            liste_sel_muscle(end+1)=sum(diff_muscle_arr(:,k) > discrim_diff);

            % If at least two Post Activation Depression per amplitude
            if sum(diff_muscle_arr(:,k) > discrim_diff) > 1
                
                first = first_pulse(:,k:k+2); % First Pulses for each paired-pulse at selected amplitude
                second = second_pulse(:,k:k+2); % Second Pulses for each paired-pulse at selected amplitude

                colormap_spec = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.4660, 0.6740, 0.1880];
                    
                figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
                
                if plot_check
                    plot (first(:,1),'--','DisplayName','First Pulse, First Repetition','Color',colormap_spec(1,:),'LineWidth',1.5)
                    hold on ;
                    plot(second(:,1),'-','DisplayName','Second Pulse, First Repetition','Color',colormap_spec(1,:),'LineWidth',1.5)
                    hold on ;
                    plot (first(:,2),'--','DisplayName','First Pulse, Second Repetition','Color',colormap_spec(2,:),'LineWidth',1.5)
                    hold on ;
                    plot(second(:,2),'-','DisplayName','Second Pulse, Second Repetition','Color',colormap_spec(2,:),'LineWidth',1.5)
                    hold on ;
                    plot (first(:,3),'--','DisplayName','First Pulse, Third Repetition','Color',colormap_spec(3,:),'LineWidth',1.5)
                    hold on ;
                    plot(second(:,3),'-','DisplayName','Second Pulse, Third Repetition','Color',colormap_spec(3,:),'LineWidth',1.5)
                    ylim([-0.5 0.5])
                    title([nom_muscles{i},' at the amplitude ',num2str(amp_values(k)),'mA'])
                    legend
                end
            end
        end
        liste_sel{end+1}=liste_sel_muscle;
    end
end


%%%%%%%%%%%


%% Calculations First and Second Normalized Activation for each Paired-Pulse

if b_rms || b_zscore
    diff_norm = [];
    for i = 1:n_muscles
       diff_muscle_norm = [];
       
       for k = 1:nb_pulse
           signal_diff_norm = EMG_cell_norm{i}(:,k);
           pulses_1_norm = signal_diff_norm(round((-t_acq_d(1)+t_stim(1))*fs):round((-t_acq_d(1)+t_stim(2))*fs));% Retrieve the signal from the first analysis window
           pulses_2_norm = signal_diff_norm(round((-t_acq_d(1)+t_stim(1)+d_int_p)*fs):round((-t_acq_d(1)+t_stim(2)+d_int_p)*fs)); % Retrieve the signal from the first analysis window
           stim_1_norm= (max(pulses_1_norm)-min(pulses_1_norm)); % Stim from first analysis window
           stim_2_norm= (max(pulses_2_norm)-min(pulses_2_norm)); % Stim from second analysis window
           diff_p_norm = stim_1_norm - stim_2_norm; % Calculation of the normalized difference for each paired-pulse (Other posibility: ratio)
           diff_muscle_norm(end+1)=diff_p_norm;
           
       end
       diff_norm(end+1,:)=diff_muscle_norm;
    end
end


%%%%%%%%%%%


%% Normalized Post Depression Activation Comparison Curve

if b_rms || b_zscore
    fig8 = figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
    
    diff_arr_norm = {};
    for i = 1:n_muscles
    
      %% Table of Normalized Difference/Ratio Calculations classified by Stimulation Amplitude
      j=1;
      a=0;
      for k = 1:length(stims)
            if range_impulsion(k) == amp_values(j) % For each Stimulation Amplitude change, the column changes
                a=a+1;
                diff_muscle_arr_norm(a,j) = diff_norm(i,k);
            else
                j=j+1;
                a=1;
                diff_muscle_arr_norm(a,j) = diff_norm(i,k);
            end
       end
       diff_arr_norm{end+1}=diff_muscle_arr_norm;
       mean_diff_norm = mean(diff_muscle_arr_norm, 'omitnan'); % Mean of normalized differences/ratios over pulses at the same amplitude
    
       key_color = strtok(nom_muscles{i}, '_');
       if i < (n_muscles/2)+1
            if plot_check
                subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
                plot(amp_values,mean_diff_norm,'--o','DisplayName',nom_muscles{i},'Linewidth',2,'Color',colors_muscles.(key_color))
                % Display Muscles Names on each curve
                text(amp_values(end) + 0.5, mean_diff_norm(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
                hold on ;
                % Plot of single points corresponding to normalized difference/ratio of each pulse
                scatter(range_impulsion,diff_norm(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
                hold on;
            end
       else
            if plot_check
                subaxis(1,2,1,'spacingHoriz',0.02,'spacingVert',0.045,'MarginRight',0.02,'MarginTop',0.03,'MarginLeft',0.03,'MarginBottom',0.07) 
                plot(amp_values,mean_diff_norm,'-o','DisplayName',nom_muscles{i},'LineWidth',2,'Color',colors_muscles.(key_color))
                % Display Muscles Names on each curve
                text(amp_values(end) + 0.5, mean_diff_norm(end), nom_muscles{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors_muscles.(key_color));
                hold on ;
                % Plot of single points corresponding to normalized difference/ratio of each pulse
                scatter(range_impulsion,diff_norm(i, 1:nb_pulse),'MarkerEdgeColor',colors_muscles.(key_color), 'HandleVisibility', 'off');
                hold on;
            end
       end
    
    end
    xlabel('Stimulation intensity (mA)')
    xlim(amp_range);
    ylim(lim_diff_norm);
    ylabel('Mean difference first and second pulse normalized amplitude  (mV)')
    if legende_m % Controls the display of legend box
        legend
    end
    title('Mean difference first and second pulse normalized amplitude over time','Color','red')
    xticks(val0:incr:amp_max)
    
    
    if save_recruit_norm
        saveas(fig8, fullfile(file_save, 'diff_norm.png'));
        saveas(fig8, fullfile(file_save, 'diff_norm.pdf'));
        saveas(fig8, fullfile(file_save, 'diff_norm.svg'));
        savefig(fig8 , fullfile(file_save, 'diff_norm.fig'));
    end
end


%%%%%%%%%%%


%% Normalized Post Activation Depression Selection Curve

if b_rms || b_zscore
    if specific_norm
        
        liste_sel_norm={};
        for i = 1:n_muscles
            
        
            first_pulse_norm = EMG_cell_norm{i}(round((t_stim(1)-t_acq_d(1))*fs):round((t_stim(2)-t_acq_d(1))*fs), 1:nb_pulse); % Contains the first pulse analysis window of each paired-pulse stimulation
            second_pulse_norm = EMG_cell_norm{i}(round((t_stim(1)+d_int_p-t_acq_d(1))*fs):round((t_stim(2)+d_int_p-t_acq_d(1))*fs), 1:nb_pulse); % Contains the second pulse analysis window of each paired-pulse stimulation
    

            diff_muscle_arr_norm = diff_arr_norm{i};
            %% Selection of paired-pulses with sufficient Normalized Post Activation Depression
            liste_sel_muscle_norm = [];
            for k = 1:length(diff_muscle_arr_norm(1,:))
                
                liste_sel_muscle_norm(end+1)=sum(diff_muscle_arr_norm(:,k) > discrim_diff_norm);
    
                % If at least two Normalized Post Activation Depression per amplitude
                if sum(diff_muscle_arr_norm(:,k) > discrim_diff_norm) > 1
         
                    first_norm = first_pulse_norm(:,k:k+2); % First Normalized Pulses for each paired-pulse at selected amplitude
                    second_norm = second_pulse_norm(:,k:k+2); % Second Normalized Pulses for each paired-pulse at selected amplitude
                        
                    figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w');
                    if plot_check
                        plot (first_norm(:,1),'--','DisplayName','First Pulse, First Repetition','Color',colormap_spec(1,:),'LineWidth',1.5)
                        hold on ;
                        plot(second_norm(:,1),'-','DisplayName','Second Pulse, First Repetition','Color',colormap_spec(1,:),'LineWidth',1.5)
                        hold on ;
                        plot (first_norm(:,2),'--','DisplayName','First Pulse, Second Repetition','Color',colormap_spec(2,:),'LineWidth',1.5)
                        hold on ;
                        plot(second_norm(:,2),'-','DisplayName','Second Pulse, Second Repetition','Color',colormap_spec(2,:),'LineWidth',1.5)
                        hold on ;
                        plot (first_norm(:,3),'--','DisplayName','First Pulse, Third Repetition','Color',colormap_spec(3,:),'LineWidth',1.5)
                        hold on ;
                        plot(second_norm(:,3),'-','DisplayName','Second Pulse, Third Repetition','Color',colormap_spec(3,:),'LineWidth',1.5)
                        ylim([-150 200])
                        title([nom_muscles{i},' at the amplitude ',num2str(amp_values(k)),'mA, Normalized'])
                        legend
                    end
                end
            end
            liste_sel_norm{end+1}=liste_sel_muscle_norm;
        end
    end
end
%%%%%%%%%%%


%% Summary Table of thresholds and ranks accross all conditions

filePath = [folderPath,'thresholds.csv'];
nb_condition = length(condition_list);

if isfile(filePath)
    resume_thresh = readcell(filePath);
else
    resume_thresh = cell(n_muscles+1, 2*nb_condition + 1); % Retrieves the information from previous conditions
    
end
resume_thresh{1, 1} = 'Nom Muscles';

% Set the headers for each condition (Line Names)
for i = 1:nb_condition
    if i <= nb_condition/2
        resume_thresh{1,2*i} = string(condition_list(i))+' Threshold';
        resume_thresh{1,2*i+1} = string(condition_list(i))+' Rank';
    else 
        resume_thresh{1,2*(i-1)+nb_condition/2} = string(condition_list(i))+' Threshold';
        resume_thresh{1,2*i-1+nb_condition/2} = string(condition_list(i))+' Rank';
    end
end

% Set the headers for each muscle (Column Names)
for i = 1:n_muscles
    resume_thresh{i+1,1} = string(nom_muscles(i));
end

% Retrieves the thresholds information for the current condition
thresh_cell = num2cell(thresh);
for i = 1:length(thresh_cell)
    if isnan(thresh_cell{i})
        thresh_cell{i} = '-'; % - if muscle never activated
    end
end

% Retrieves the normalized thresholds information for the current condition
thresh_norm_cell = num2cell(thresh_norm);
for i = 1:length(thresh_norm_cell)
    if ~(b_rms || b_zscore)
        thresh_norm_cell{i} = '-'; % - if muscle never activated
    else
        if isnan(thresh_norm_cell{i})
            thresh_norm_cell{i} = '-'; % - if muscle never activated
        end
    end
end

% Retrieves the ranks information for the current condition
rank_cell = num2cell(rank);
for i = 1:length(rank_cell)
    if isnan(rank_cell{i}) 
        rank_cell{i} = '-'; % - if muscle never activated
    end
end

% Retrieves the normalized ranks information for the current condition
rank_norm_cell = num2cell(rank_norm);
for i = 1:length(rank_norm_cell)
    if ~(b_rms || b_zscore)
        rank_norm_cell{i} = '-'; % - if muscle never activated
    else
        if isnan(rank_norm_cell{i})
            rank_norm_cell{i} = '-'; % - if muscle never activated
        end
    end
end

% Adds current condition calculations to the table
for i = 2:n_muscles+1
    resume_thresh{i,2*num_condition} = thresh_cell{i-1};
    resume_thresh{i, nb_condition + 2*num_condition} = thresh_norm_cell{i-1};
    resume_thresh{i,1+ nb_condition/2 + num_condition} = rank_cell{i-1};
    resume_thresh{i,1+ num_condition + nb_condition*(3/2)} = rank_norm_cell{i-1};
end

writecell(resume_thresh,filePath)


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


%%%%%%%%%%%


%% Display Muscle Button

function state = boolToOnOff(val)
    if val
        state = 'on';
    else
        state = 'off';
    end
end




