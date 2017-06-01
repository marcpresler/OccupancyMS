% OccupancyCode_Superscript_v1.m
% Marc Presler, Martin Wuehr, Allon Klein, Elizabeth Van Itallie 
% December 16th 2016


% Brief explanation of approach: 
%   Stoichiometry of a phospho-site can be represented as a slope which is
%   set by the corresponding changes in the phos and unmodified form of a
%   species. From the input data, the algorithm fits this slope to calculate 
%   site stoichiometry. Bootstrapping is performed to establish confidence intervals.

% Script returns phospho-site occupancy from the input of the matching phospho 
%   and unmodified trends from quantitative mass spectrometry data. User must
%   make sure that the information in the Data Import section is correct for
%   their data. 

% Table of contents:
% 1) Data import  
% 2) Parameters
% 3) Estimate Site Occupancy and Determine CIs with Bootstrapping
% 4) Plot Data 
% 5) Consolidate and export data




%% 1) Data import  

% Loads: 
% .xlxs file containing TMT signal for phos sites and matched unmodified peptides
% Mean normalization of the data is not necesssary, but helpful for
% plotting.

% Column names in .xlsx file correpsonding to unmodified and phos-sites data are input by the users 
% in 'peptide_data_variables' and 'phospho_data_variables'.
% '.xlsx' file requires a column with 'GeneSymbol' and 'SitePosition'.

    % Enter filename for import as a string with file extension, e.g., 'SampleFile_v1.xlsx'
filename = 'SampleFile_v1.xlsx';

    % Enter variable names of the Non-phospho peptide data (cell array)
peptide_data_variables = {'rq_126_sn','rq_127n_sn','rq_127c_sn','rq_128n_sn','rq_128c_sn','rq_129n_sn','rq_129c_sn','rq_130n_sn','rq_130c_sn','rq_131_sn'};

    % Enter variable names of the phospho peptide data (cell array)
phospho_data_variables = {'rq_126_Phos_Normalized_Mean_sn','rq_127n_Phos_Normalized_Mean_sn','rq_127c_Phos_Normalized_Mean_sn','rq_128n_Phos_Normalized_Mean_sn','rq_128c_Phos_Normalized_Mean_sn','rq_129n_Phos_Normalized_Mean_sn','rq_129c_Phos_Normalized_Mean_sn','rq_130n_Phos_Normalized_Mean_sn','rq_130c_Phos_Normalized_Mean_sn','rq_131_Phos_Normalized_Mean_sn'};

    % Enter additional labeling information. Code requires 'GeneSymbol' and 'SitePosition'
labeling_variables =  {'ProteinId','GeneSymbol','Description','SitePosition', 'site_id', 'sequence'};

    % Determines number of conditions
what_plex_is_data = size(peptide_data_variables,2);


%Imports data from spreadsheet into the form of a table 
single_Phos_and_NON_phos = readtable(filename); 
    length_of_dataset = size(single_Phos_and_NON_phos,1);


%% 2) Parameters

% Confidence interval Parameters: 

    % Set the target confidence interval between 0 and 1 (e.g., for 95% CI, enter 0.95) 
target_ci = 0.95;
    
    % Set number of times to bootstrap. 1000-10,000 recommend, but code becomes
    % quite slow for more than a handful of sites if iterations > 100. 
bootci_iterations = 1000;


%Other options:

    % Enter 1 to save .xlsx file of results. Filename with parameters is generated automatically  
save_table = 1;

    % Enter 1 to see the fitting plots used to calculate stoichiometry. If
    % dataset exceeds ~10 sites, set to zero as Matlab may crash due
    % to a graphics error. 
plot_fitting_scatter = 1;

    % Enter 1 to plot occupancy trends. Set x-axis and color option below.
    % Figure 1 is the unmodified versus Phos Trends
    % Figure 2 is the subplot of occupancy trends
    % Figure 3 is occupancy trends plotted together.
    % With large dataset, consider setting to zero 
plot_data = 1;

    % Enter appropriate x axis for plots
    x_axis_data = 0:2:18;
    
    % Enter color
    colors = 'c';

    %Set mean width of confidence interval in percent. Set to 100 to include all data 
 CI_width_cutoff = 100;

    % "Corrected" sites are a small subset of the data that are high quality and 
    % are clearly low or high occupancy trends, but give unrealistically wide
    % confidence intervals. This can occur because small amounts of error as 
    % slopes approach zero or infinity (i.e., towards 0 or toward 100% occupancy)
    % can dominate to give extreme error during bootstrapping. 
    % For these sites, we set the CI's arbitarily at 0 to 25%, or 75 to 100%. 
    % Type 1 to remove these from the dataset.  
    removed_all_corrected_sites = 0;


%Sets consistant random number set
rng('default')


%% 3) Estimate Site Occupancy and Determine CIs with Bootstrapping

%preallocate variables 
        Fit_confidence_lower_stored_single = zeros(size(single_Phos_and_NON_phos,1),size(phospho_data_variables,2));
        occupancytrend_single = zeros(size(single_Phos_and_NON_phos,1),size(phospho_data_variables,2));
        Fit_confidence_higher_stored_single = zeros(size(single_Phos_and_NON_phos,1),size(phospho_data_variables,2));

        section_2_indicies = zeros(length_of_dataset,1);
        section_1_indicies = zeros(length_of_dataset,1);
        section_3_indicies = zeros(length_of_dataset,1);

%For each sites, calculate the stoichiometry and confidence intervals 
for site_counter = 1:length_of_dataset
    
                %Used to make fitting figures, if set to 1 above
                if plot_fitting_scatter
                    figure
                end
                
   
    %Store raw data for one site per loops
    Phos_Single = single_Phos_and_NON_phos{site_counter,phospho_data_variables};
    NON_phos = single_Phos_and_NON_phos{site_counter,peptide_data_variables};
    
    %Code will fail if there are zeros in channels, since it may divide by zero later on. 
    %Replace zeros with 1E-9 for both phos and nonphos
    if sum(Phos_Single==0)>=1
       Phos_Single(logical(Phos_Single==0))=1E-9;
    end
    
    if sum(NON_phos==0)>=1
       NON_phos(logical(NON_phos==0))=1E-9;
    end
    
    
    % If a trend contains identical points (e.g., multiple zeros) in one species, a
    % subsample  will at some frequency give an infinite or 0 slope and fail during bootstrapping. 
    % This flag allows the downstream function to avoid this error. 
    FLAG_repeated_points = 0;
    if size(unique(Phos_Single),2) < what_plex_is_data || size(unique(NON_phos),2) < what_plex_is_data
        FLAG_repeated_points = 1;
    end
    
    
   
    %Create matrix of normalized to each condition...Will use this to plot the 10 ratios for each and plot them against each other 
        %For single phos
    [columnrepmat_Phos_Single,rowrepmat_Phos_Single] = meshgrid(Phos_Single, Phos_Single);
    ratios_Phos_Single = (columnrepmat_Phos_Single./rowrepmat_Phos_Single);
        %for Non-phos 
    [columnrepmat_NON_phos,rowrepmat_NON_phos] = meshgrid(NON_phos, NON_phos);
    ratios_NON_phos = (columnrepmat_NON_phos./rowrepmat_NON_phos);


    %Estimate phos/nonphos ratios with fitting
            %preallocate
        Fit_occupancy_single_optimal_stored = zeros(1,what_plex_is_data);
        Fit_confidence_lower_stored = zeros(1,what_plex_is_data);
        Fit_confidence_higher_stored = zeros(1,what_plex_is_data);
        
        %Calculate optimal solution  
        for condition_counter_optimal = 1:what_plex_is_data
            %call the corresponding non phos and phos reference point normalized for optimization
            refPoint_Normalized_NONphosOptimal = ratios_NON_phos(condition_counter_optimal,:);
            refPoint_Normalized_PhosOptimal = ratios_Phos_Single(condition_counter_optimal,:);
            %reshape matrix to 10,2 for compatablility with bootci input format
            input_matrix=[refPoint_Normalized_NONphosOptimal;refPoint_Normalized_PhosOptimal]';
            
            %Perform fitting, orthogonal transformation, calculates percent occupancy, and corrects for 
            %"impossible" values above or below zero %
            [Fit_occupancy_optimal] = call_TLS_fitting_v1(input_matrix, FLAG_repeated_points);
            
            %Performs the same function as above but with no correction.
            %This is useful for evaluating small number of "correctable" 
            %sites which are clearly qualitatively low or high, but the 
            %answer gives "borderline" impossible values due to noise. In these cases, the answers will 
            %be just barely below zero or barely above 100 without correction.
            %These cases are set to 0 or 100, with CI's arbitrarily set to 0 to 25, or 75 to 100%. 
            Fit_occupancy_optimal_NoCorrection = call_TLS_fitting_NO_corrections_v1(input_matrix);
            
                Fit_occupancy_single_optimal_stored(:,condition_counter_optimal) = Fit_occupancy_optimal(1,1);

                Fit_occupancy_non_optimal_stored(:,condition_counter_optimal) = Fit_occupancy_optimal(2,1);
                
                Fit_occupancy_single_optimal_NOcorrection_stored(:,condition_counter_optimal) = Fit_occupancy_optimal_NoCorrection(1,1);

        end
                
        
        %Calculate confidence interval by bootstrapping. Evaluates if optimal estimate is within 0 to 100 bounds, where bootstrapping will proceed.   
        if max(Fit_occupancy_single_optimal_NOcorrection_stored) > 0 && min(Fit_occupancy_single_optimal_NOcorrection_stored) < 100
                        section_boostrap_num(site_counter) = 1;

            for timepoint_counter_CIs = 1:what_plex_is_data
                %again call the corresponding non phos and phos reference point normalized for optimization
                refPoint_Normalized_NONphos = ratios_NON_phos(timepoint_counter_CIs,:);
                refPoint_Normalized_Phos = ratios_Phos_Single(timepoint_counter_CIs,:);
                %reshape matrix to 10,2 for compatablility with bootci input format
                input_matrix=[refPoint_Normalized_NONphos;refPoint_Normalized_Phos]';

                %reshape matrix to 10,2 for compatablility with bootci input format
                input_matrix=[refPoint_Normalized_NONphos;refPoint_Normalized_Phos]';
               
                ci_occupancies = bootci(bootci_iterations,{@call_TLS_fitting_v1,input_matrix,FLAG_repeated_points},'alpha',(1-target_ci),'type','bca');


                
                
               
                    %Plots the line fitting for each othe timepoints used
                    %to calculate stoichiometry. Optimal result is plotted
                    %in red, the confidence intervals are in blue. The 
                    %reference (or the point the line is normalized to) is
                    %always (1,1), which is marked in black. The axes
                    %may require adjustment for best visualization
                if plot_fitting_scatter
                   
                    %calc optimal vector for the purposes of plotting...
                    % P gives coefficiencts used to plot
                    % Yhat1 gives the optimal vector from TLS regression 
                    [~, P, Yhat1] = fit_2D_data_modified_v1(refPoint_Normalized_NONphos,refPoint_Normalized_Phos,'no');

                    ratio_low = -1./((ci_occupancies(1,1)/100)/(1-(ci_occupancies(1,1)/100)));
                    ratio_high = -1./((ci_occupancies(2,1)/100)/(1-(ci_occupancies(2,1)/100)));
                    subplot(4,3,timepoint_counter_CIs)
                    plot(refPoint_Normalized_NONphos,refPoint_Normalized_Phos,'.','color',[0.5 0.5 0.5],'markersize',15); 
                    hold on;
                    plot(refPoint_Normalized_NONphos,(Yhat1),'r-');
                    hold on
                    plot(1,1,'.','Markersize',20,'color','k')
                    hold on
                    plot(refPoint_Normalized_NONphos, refPoint_Normalized_NONphos*ratio_low + P(1,2),'c-');
                    hold on
                    plot(refPoint_Normalized_NONphos, refPoint_Normalized_NONphos*ratio_high + P(1,2),'c-');
                    set(gca,'fontsize',12)
                    title(['Ref. Point ',num2str(timepoint_counter_CIs)],'fontsize',10)
                    xlabel('NonPhos Ratios','fontsize',10)
                    ylabel('Phos Ratios','fontsize',10)
                    axis square 
                end

                
                
                
                
                %Store lower confidence for full time series, phospho
                Fit_confidence_lower_single(:,timepoint_counter_CIs) = ci_occupancies(1,1);

                %Store higher confidence, for full time series, phospho
                Fit_confidence_higher_single(:,timepoint_counter_CIs) = ci_occupancies(2,1);

                %Store lower confidence, for full time series, nonphos
                Fit_confidence_lower_non(:,timepoint_counter_CIs) = ci_occupancies(1,2);
              
                %Store lower confidence, for full time series, nonphos
                Fit_confidence_higher_non(:,timepoint_counter_CIs) = ci_occupancies(2,2);

            end
            
        %assess if in section 3, which are sites with no certainty due to excessive noise. Sets CIs to 0 and 100%. Bootrapping not necessary here.   
        elseif max(Fit_occupancy_single_optimal_NOcorrection_stored) < -2  || max(Fit_occupancy_single_optimal_NOcorrection_stored) > 102
            section_3_indicies(site_counter,1) = 1;
            
            Fit_confidence_lower_single = zeros(1,what_plex_is_data);
            
            Fit_confidence_higher_single = repmat(100,1,what_plex_is_data);
            
            Fit_confidence_lower_non = zeros(1,what_plex_is_data);
            Fit_confidence_higher_non = repmat(100,1,what_plex_is_data);
            
        %"Correctable Sites", assess if in "section 1," which are sites that are clearly low occupancy qualitatively but still give marginally nonsense values (i.e., -1% occupancy). 
        % Set optimal to all 0, and CIs to 0 to 25%. 
        % Can remove these with option at beginning of code. This will be a small percentage of the data.
        elseif max(Fit_occupancy_single_optimal_NOcorrection_stored) > -2  && max(Fit_occupancy_single_optimal_NOcorrection_stored) < 0
            section_1_indicies(site_counter,1) = 1;            
            
            Fit_confidence_lower_single = zeros(1,what_plex_is_data);
            
            Fit_confidence_higher_single = repmat(25,1,what_plex_is_data);
            
            Fit_confidence_lower_non = zeros(1,what_plex_is_data);
            Fit_confidence_higher_non = repmat(100,1,what_plex_is_data);
            
          %assess if in section 2, which are, sites that are clearly 100% but give large errors as above. In practice, this is less common than class 1.    
        elseif max(Fit_occupancy_single_optimal_NOcorrection_stored) > 100  && max(Fit_occupancy_single_optimal_NOcorrection_stored) < 101
            section_2_indicies(site_counter,1) = 1;
            
            Fit_confidence_lower_single = repmat(75,1,what_plex_is_data);
            
            Fit_confidence_higher_single = repmat(100,1,what_plex_is_data);
            
            Fit_confidence_lower_non = zeros(1,what_plex_is_data);
            Fit_confidence_higher_non = repmat(100,1,what_plex_is_data);      
        

        end
            %Store all the occupancies and their confidence bounds for each site

                Fit_confidence_lower_stored_single(site_counter,:) = Fit_confidence_lower_single;

                Fit_confidence_lower_stored_non(site_counter,:) = Fit_confidence_lower_non;


                occupancytrend_single(site_counter,:)=Fit_occupancy_single_optimal_stored;

                occupancytrend_non(site_counter,:)=Fit_occupancy_non_optimal_stored;


                Fit_confidence_higher_stored_single(site_counter,:) = Fit_confidence_higher_single;

                Fit_confidence_higher_stored_non(site_counter,:) = Fit_confidence_higher_non;
          
end

%accounts for all corrected sites
indicies_all_corrected_sites = section_2_indicies + section_1_indicies + section_3_indicies;
num_Corrected_sites = sum(indicies_all_corrected_sites);
num_section_3_indicies = sum(section_3_indicies);
num_section_1_indicies = sum(section_1_indicies);
num_section_2_indicies = sum(section_2_indicies);
%accounts for sites bootrapped (vast majority)
section_boostrap_num = sum(section_boostrap_num);
    
     
     
%% 4) Plot Data 

% If data is to be plotted, input data must be normalized (e.g., to the trend mean) 
% before running the script.

    % Enter 1 to plot occupancy trends. Set x-axis and color option below.
    % Figure 1 is the unmodified versus Phos Trends
    % Figure 2 is the subplot of occupancy trends
    % Figure 3 is occupancy trends plotted together.
    % With large dataset, set to zero. 

x_axis_labeling = 'Time (min)';

if plot_data  
        length_of_dataset = size(single_Phos_and_NON_phos,1);
        dim_subplot_first = round(sqrt(length_of_dataset));
        dim_subplot_second = round(length_of_dataset/dim_subplot_first);
        if (dim_subplot_first*dim_subplot_second)<length_of_dataset
            dim_subplot_second = dim_subplot_second + 1;
        end
       figure
       for site_counter = 1:length_of_dataset;
        subplot(dim_subplot_first,dim_subplot_second,site_counter)
        plot(x_axis_data,single_Phos_and_NON_phos{site_counter,phospho_data_variables},'-','color','g','LineWidth',2,'marker','.','markersize',12)
        hold on
        plot(x_axis_data,single_Phos_and_NON_phos{site_counter,peptide_data_variables},'-','color','b','LineWidth',2,'marker','.','markersize',12)
        hold off
        box off
        xlabel(x_axis_labeling)
        ylabel('Relative Abundance')
        sitenum = single_Phos_and_NON_phos{site_counter,'SitePosition'};
        pretitle = single_Phos_and_NON_phos{site_counter,'GeneSymbol'};
        title(pretitle,'fontsize',10)
        set(gca,'fontsize',10)
        legend(['Phospho-Site ',num2str(sitenum)],'NonPhos')
       end
       set(gcf,'color','w')
       
     
     
    %with shaded plots and not plotting pptase...
           figure
        dim_subplot_first2 = round(sqrt(length_of_dataset));
        dim_subplot_second2 = round(length_of_dataset/dim_subplot_first2);
        if (dim_subplot_first2*dim_subplot_second2)<length_of_dataset
            dim_subplot_second2 = dim_subplot_second2 + 1;
        end

     for site_counter = 1:length_of_dataset;
        x_axis_data = 0:2:18;
        upper = Fit_confidence_higher_stored_single(site_counter,:);
        lower = Fit_confidence_lower_stored_single(site_counter,:);
        optimal_trend = occupancytrend_single(site_counter,:);
        
        subplot(dim_subplot_first2,dim_subplot_second2,site_counter)

        fill([x_axis_data, fliplr(x_axis_data)], [upper, fliplr(lower)], colors, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on      
        plot(x_axis_data, optimal_trend,'color',colors,'linestyle','-', 'linewidth',1.5)
        hold on
         
        sitenum = single_Phos_and_NON_phos{site_counter,'SitePosition'};
        pretitle2 = single_Phos_and_NON_phos{site_counter,'GeneSymbol'};
%         title(strcat(pretitle2,' Site ',sitenum),'fontsize',18)    
        title([pretitle2{1,1},', Site ',num2str(sitenum)],'fontsize',12) 
        xlabel(x_axis_labeling)
        ylabel('Percent Occupancy')
        ylim([0,100])
        set(gca,'fontsize',14)
        box off
     end
     set(gca,'Layer', 'Top')
     set(gcf,'color','w')

     
     
         %with shaded plots and not plotting pptase multiplot...
           figure
      for site_counter = 1:length_of_dataset;
        x_axis_data = 0:2:18;
        upper = Fit_confidence_higher_stored_single(site_counter,:);
        lower = Fit_confidence_lower_stored_single(site_counter,:);
        optimal_trend = occupancytrend_single(site_counter,:);
        
        fill([x_axis_data, fliplr(x_axis_data)], [upper, fliplr(lower)], colors, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on      
        plot(x_axis_data, optimal_trend,'color',colors,'linestyle','-', 'linewidth',1.5)
        hold on
         
        sitepositionnum = single_Phos_and_NON_phos{site_counter,'SitePosition'};
        pretitle2 = single_Phos_and_NON_phos{site_counter,'GeneSymbol'};
        title('Site occupancy over time','fontsize',12)    
        text(x_axis_data(end)+0.25,optimal_trend(end)+2,[pretitle2{1,1},', Site ',num2str(sitenum)],'fontsize',12,'color','k')     
        xlabel(x_axis_labeling)
        ylabel('Percent Occupancy')
        ylim([0,100])
        set(gca,'fontsize',14)
        box off
      end
      set(gca,'Layer', 'Top')
    x_axis_data = 0:2:18;
    ax = gca;
    c = ax.PlotBoxAspectRatio;
    ax.PlotBoxAspectRatio = [0.7043    1.0000    0.7043];
    set(gcf,'color','w')


end 

%% 5) Consolidate and export data

% Returns:      
%   The code exports a table of  occupancies with the upper and lower bound 
% of the confidence intervals at a given cutoff, if specified. 

%make stored dataset of optimal occupancy estimates
    %preallocate
    occupancytrend_array = zeros(length_of_dataset,what_plex_is_data);
for site_counter = 1:length_of_dataset; 
    occupancytrend_array(site_counter,:) = occupancytrend_single(site_counter,:);
    occupancytrend_single_labeling(site_counter,:) = single_Phos_and_NON_phos(site_counter,labeling_variables);
end


%make stored dataset of higher confidence bound
    %preallocate
    Fit_confidence_higher_stored_single_array = zeros(length_of_dataset,what_plex_is_data);
for site_counter = 1:length_of_dataset; 
    Fit_confidence_higher_stored_single_array(site_counter,:) = Fit_confidence_higher_stored_single(site_counter,:);
end

%make stored dataset of lower confidence bound 
    %preallocate
    Fit_confidence_lower_stored_single_array = zeros(length_of_dataset,what_plex_is_data);
for site_counter = 1:length_of_dataset; 
    Fit_confidence_lower_stored_single_array(site_counter,:) = Fit_confidence_lower_stored_single(site_counter,:);
end


% Calculates width of CI at each time point for every trend 
CI_distance = minus(Fit_confidence_higher_stored_single,Fit_confidence_lower_stored_single);

% Index of sites which pass threshold for confidence 
pass_cuttoff_indices = mean(CI_distance,2)<CI_width_cutoff;


% Removes "corrected sites" from dataset if set at 1 above
if removed_all_corrected_sites && bootci_iterations > 0
    indicies_Uncorrected_Sites = indicies_all_corrected_sites==0;

    pass_cuttoff_indices = (pass_cuttoff_indices + indicies_Uncorrected_Sites)==2;
end


%Converts array of data into tables for each trend and the CIs 
occupancytrend_array = [occupancytrend_array,mean(occupancytrend_array,2)];
occupancytrend_array_table = array2table(occupancytrend_array);
occupancytrend_array_table.Properties.VariableNames = {'occ_126_phospho','occ_127n_phospho','occ_127c_phospho','occ_128n_phospho','occ_128c_phospho','occ_129n_phospho','occ_129c_phospho','occ_130n_phospho','occ_130c_phospho','occ_131_phospho','mean_occupancy'};

Fit_confidence_higher_stored_single_array_table = array2table(Fit_confidence_higher_stored_single_array);
Fit_confidence_higher_stored_single_array_table.Properties.VariableNames = {'highCI_126_phospho','highCI_127n_phospho','highCI_127c_phospho','highCI_128n_phospho','highCI_128c_phospho','highCI_129n_phospho','highCI_129c_phospho','highCI_130n_phospho','highCI_130c_phospho','highCI_131_phospho'};

Fit_confidence_lower_stored_single_array_table = array2table(Fit_confidence_lower_stored_single_array);
Fit_confidence_lower_stored_single_array_table.Properties.VariableNames = {'lowCI_126_phospho','lowCI_127n_phospho','lowCI_127c_phospho','lowCI_128n_phospho','lowCI_128c_phospho','lowCI_129n_phospho','lowCI_129c_phospho','lowCI_130n_phospho','lowCI_130c_phospho','lowCI_131_phospho'};


%Adds relative trends to be exported too
A_optimal_occupancy_table = [occupancytrend_single_labeling,occupancytrend_array_table...
    Fit_confidence_higher_stored_single_array_table...
    Fit_confidence_lower_stored_single_array_table...
    single_Phos_and_NON_phos(:,phospho_data_variables)...
    single_Phos_and_NON_phos(:,peptide_data_variables)];

A_optimal_occupancy_table_cutoff = A_optimal_occupancy_table(pass_cuttoff_indices,:);


%Automated labeling of filename with bootstrapping parameters
bootstrap_num_label = num2str(bootci_iterations);
date_format = 'yyyymmdd';
date_label = datestr(now,date_format);
stoich_table_filename = [date_label,'_Occupancy_trends_','ci_',num2str(target_ci*100),'_',bootstrap_num_label,'x_from_',filename];

%save table
if save_table
    writetable(A_optimal_occupancy_table_cutoff,stoich_table_filename);
end

