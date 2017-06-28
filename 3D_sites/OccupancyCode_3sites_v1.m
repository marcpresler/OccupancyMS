% OccupancyCode_3sites_v1.m
% Marc Presler, Elizabeth Van Itallie, Martin Wuehr, Allon Klein 
% June 28rd, 2017


% Brief explanation of approach: 
%   Stoichiometry of a phospho-site can be represented as a slope which is
%   set by the corresponding changes in the phos and unmodified form of a
%   species. From the input data, the algorithm fits this slope to calculate 
%   site stoichiometry in multidimensional space, implemented here in 3D using
%   Singular Value Decomposition/Principal Component Analysis. 
%   Bootstrapping is performed to establish confidence intervals.

% Script returns phospho-site occupancy from the input of the matching phospho 
%   and unmodified trends from quantitative mass spectrometry data. User must
%   modify the Set Data section for their own data.  

% Table of contents:
% 1) Set Data
% 2) Parameters
% 3) Estimate Site Occupancy and Determine CIs with Bootstrapping
% 4) Plot Data 


%% Set Data

%Example dataset with good confidence.
protein_name = 'CamII kinase gamma subunit';
Unmodified = [1.853826761	1.210046793	0.772253604	0.664231207	0.641812548	0.926143842	0.77349089	0.940596628	1.060760367	1.156837361];
Phos_form1 = [0.247688829	0.819026015	1.126183959	1.271421608	1.313246744	1.207645023	1.00520424	1.063844724	0.952795287	0.992943572];
Phos_form2 = [0.078187751	1.16440218	1.723378104	1.639999719	1.340729896	1.060284057	0.917790465	0.820619092	0.668347365	0.586261371];

%Example with poor confidence (Phos forms are too similar). Uncomment to use. 
% protein_name = 'MAP4';
% Unmodified = [0.642625763	0.824391926	0.73422287	0.783439989	0.932059961	0.960156128	1.053610222	1.170085736	1.44752762	1.451879785];
% Phos_form1 = [1.375295433	1.241275277	1.252750086	1.131360476	0.976467189	0.932912503	0.970358453	0.824224074	0.722300557	0.573055953];
% Phos_form2 = [1.316880527	1.138411643	1.219133455	1.112183659	0.987982509	0.948369948	0.924660611	0.848880648	0.784433784	0.719063216];

%% 2) Parameters

% Set the target confidence interval between 0 and 1 (e.g., for 90% CI, enter 0.90) 
target_ci = 0.90;

% Set number of times to bootstrap. 1000-10,000 recommend.
boot_iterations = 1000;

%Other settings and options:

% Enter 1 to plot occupancy trends. Set x-axis option below.
plot_occupancy = 1;

    %set appropriate x_axis for plotting, should match number of conditions
    x_axis_values = 0:2:18;
    
% Enter 1 to plot visualize plan fittings.  
plot_fitting_scatter = 1;

    %Ensures proper figure generation.   
    if plot_fitting_scatter                 
        figure
    end

%Sets consistant random number set
rng('default')

%% 3) Estimate Site Occupancy and Determine CIs with Bootstrapping

%Create all ratios, normalized to each timepoint respectively...Will use this to plot the 10 ratios for each and plot them against each other 
    %For single phos
    [columnrepmat_Phos_Single,rowrepmat_Phos_Single] = meshgrid(Phos_form1, Phos_form1);
    ratios_Phos_Single = (columnrepmat_Phos_Single./rowrepmat_Phos_Single);
    %for Non-phos 
    [columnrepmat_unmodified,rowrepmat_unmodified] = meshgrid(Unmodified, Unmodified);
    ratios_unmodified = (columnrepmat_unmodified./rowrepmat_unmodified);

    %For first phos form
    [columnrepmat_Phos_Double,rowrepmat_Phos_Double] = meshgrid(Phos_form2, Phos_form2);
    ratios_Phos_Double = (columnrepmat_Phos_Double./rowrepmat_Phos_Double);

        %Preallocate variables 
        Optimal_occupancy_trend_Stored = zeros(3,size(x_axis_values,2));
        
        occupancy_lower_stored_Unmodified = zeros(1,size(x_axis_values,2));
        occupancy_higher_stored_Unmodified = zeros(1,size(x_axis_values,2));
        
        occupancy_lower_stored_phosform1 = zeros(1,size(x_axis_values,2));
        occupancy_higher_stored_phosform1 = zeros(1,size(x_axis_values,2));
        
        occupancy_lower_stored_phosform2 = zeros(1,size(x_axis_values,2));   
        occupancy_higher_stored_phosform2 = zeros(1,size(x_axis_values,2));     
         
%Estimate Occupancy per timepoint using Singular Value Decomposition to
%estaimte the orthogonally fit plane. 
        for timepoint_counter = 1:size(x_axis_values,2)
            
            %Create matrix of normalized to each condition...Will use this to plot the 10 ratios for each and plot them against each other 
            refPoint_Normalized_NONphos = ratios_unmodified(timepoint_counter,:);
            refPoint_Normalized_Phos = ratios_Phos_Single(timepoint_counter,:);
            refPoint_Normalized_Double_Phos = ratios_Phos_Double(timepoint_counter,:);
            
            input_matrix_3D = [refPoint_Normalized_NONphos;refPoint_Normalized_Phos;refPoint_Normalized_Double_Phos]';
            
            %Estimate occupany of each form using SVD 
            optimal_occupancy_trend = call_fit_3D_data_SVD_v1(input_matrix_3D);
            
            %Obtain confidence interval through bootstrapping
            ci_occupancies = bootci(boot_iterations,{@call_fit_3D_data_SVD_v1,input_matrix_3D},'alpha',(1-target_ci),'type','bca');
            
                    %Visualizing the plane defined by the first two principal
                    %components (red) and the orthogonal vector (blue), which contains the solution. 
                    if plot_fitting_scatter                 
                      subplot(4,3,timepoint_counter)
                        [normal_vector, MeanCenteredData,Vmatrix_3d, Smatrix_3d] = fit_3D_data_modified_SVD_v1(refPoint_Normalized_NONphos',refPoint_Normalized_Phos',refPoint_Normalized_Double_Phos');

                        plot3(MeanCenteredData(:,1),MeanCenteredData(:,2),MeanCenteredData(:,3),'.','markersize',14,'color',[0.6 0.6 0.6]);
                        hold on

                        range = -2:1:2;
                        plot3(normal_vector(1)*range,normal_vector(2)*range,normal_vector(3)*range,'--','color','c','linewidth',1)

                        PC1_vector = quiver3(0,0,0,Vmatrix_3d(1,1),Vmatrix_3d(2,1),Vmatrix_3d(3,1),Smatrix_3d(1,1),'linewidth',1.5,'color','r');
                        PC1_vector_minus = quiver3(0,0,0,-Vmatrix_3d(1,1),-Vmatrix_3d(2,1),-Vmatrix_3d(3,1),Smatrix_3d(1,1),'linewidth',1.5,'color','r');

                        PC2_vector = quiver3(0,0,0,Vmatrix_3d(1,2),Vmatrix_3d(2,2),Vmatrix_3d(3,2),Smatrix_3d(2,2),'linewidth',1.5,'color','r');
                        PC2_vector_minus = quiver3(0,0,0,-Vmatrix_3d(1,2),-Vmatrix_3d(2,2),-Vmatrix_3d(3,2),Smatrix_3d(2,2),'linewidth',1.5,'color','r');

                        PC3_vector = quiver3(0,0,0,Vmatrix_3d(1,3),Vmatrix_3d(2,3),Vmatrix_3d(3,3),Smatrix_3d(3,3),'linewidth',1.5,'color','c');
                        PC3_vector_minus = quiver3(0,0,0,-Vmatrix_3d(1,3),-Vmatrix_3d(2,3),-Vmatrix_3d(3,3),Smatrix_3d(3,3),'linewidth',1.5,'color','c');

                        PC1_values = [PC1_vector.UData, PC1_vector.VData, PC1_vector.WData];
                        PC2_values = [PC2_vector.UData, PC2_vector.VData, PC2_vector.WData];

                        X_grid = [PC1_values(1).*Smatrix_3d(1,1), PC2_values(1).*Smatrix_3d(2,2),;...
                            -PC2_values(1).*Smatrix_3d(2,2), -PC1_values(1).*Smatrix_3d(1,1)];

                        Y_grid = [PC1_values(2).*Smatrix_3d(1,1), PC2_values(2).*Smatrix_3d(2,2),;...
                            -PC2_values(2).*Smatrix_3d(2,2), -PC1_values(2).*Smatrix_3d(1,1)];

                        Z_grid = (-(normal_vector(1).*X_grid + normal_vector(2).*Y_grid))./(normal_vector(3));
                        mesh(X_grid,Y_grid,Z_grid,...
                            'EdgeAlpha',0,...
                            'FaceColor','r',...
                            'FaceAlpha',0.5)
                        hold off
                        grid on
%                         campos([24 -24 10])
                        campos([15 -4.5 3])
                        set(gca,'fontsize',10)
                        title(['Ref. Point ',num2str(x_axis_values(timepoint_counter)),'min'],'fontsize',10)               
                        xlabel('Unmodified','fontsize',10)
                        ylabel('Phos','fontsize',10)
                        zlabel('Double Phos','fontsize',10)
                        
                        h = get(gca,'DataAspectRatio'); 
                        set(gca,'DataAspectRatio',[1 1 1])

                        upperaxis = max(max([X_grid; Y_grid;Z_grid]));
                        loweraxis = -upperaxis;
                        xlim([loweraxis, upperaxis])
                        ylim([loweraxis, upperaxis])
                        zlim([loweraxis, upperaxis])
                    end       
            
            %Store data per each time point
            Optimal_occupancy_trend_Stored(:,timepoint_counter) = optimal_occupancy_trend;
                        
            occupancy_lower_stored_Unmodified(:,timepoint_counter) = ci_occupancies(1,1);
            occupancy_higher_stored_Unmodified(:,timepoint_counter) = ci_occupancies(2,1);
            
            occupancy_lower_stored_phosform1(:,timepoint_counter) = ci_occupancies(1,2);
            occupancy_higher_stored_phosform1(:,timepoint_counter) = ci_occupancies(2,2);
            
            occupancy_lower_stored_phosform2(:,timepoint_counter) = ci_occupancies(1,3);   
            occupancy_higher_stored_phosform2(:,timepoint_counter) = ci_occupancies(2,3);
            

            
        end
       
        %Pull out occupancy of each form
        occupancytrend_unmodified = Optimal_occupancy_trend_Stored(1,:);
        occupancytrend_single = Optimal_occupancy_trend_Stored(2,:);
        occupancytrend_double = Optimal_occupancy_trend_Stored(3,:);
            
       
%% 4) Plot Data 

if plot_occupancy
    
    markersizeforplots = 20;
    graphfontsize = 10;
    figure 
    subplot(1,4,1)
    plot(x_axis_values,Phos_form1,'g-','LineWidth',2,'marker','.','markersize',markersizeforplots)
    hold on
    plot(x_axis_values,Unmodified,'b-','LineWidth',2,'marker','.','markersize',markersizeforplots)
    hold on
    plot(x_axis_values,Phos_form2,'color',[1 .5 0],'LineWidth',2,'marker','.','markersize',markersizeforplots)
    hold off
    box off
    xlabel('Time (min)')
    ylabel('Relative Abundance','fontsize',15)
    text(x_axis_values(end)+0.5,Phos_form1(end),'P^1','fontsize',18,'color','g')
    text(x_axis_values(end)+0.5,Phos_form2(end),'P^2','fontsize',18,'color',[1 .5 0])
    text(x_axis_values(end)-2,Unmodified(end)+0.1,'U','fontsize',18,'color','b')
    set(gca,'fontsize',graphfontsize)
    
    subplot(1,4,2)
        plot(x_axis_values,occupancy_lower_stored_Unmodified,'w','linewidth',0.5)
  
        h2 = fill([0, x_axis_values, x_axis_values(end)],...        
          [100, occupancy_lower_stored_Unmodified, 100],...
          'b','EdgeColor','none');
        set(h2,'facealpha',.15)
        hold on
        plot(x_axis_values, occupancytrend_unmodified,'b-', 'linewidth',1)%,'marker','.','markersize',markersizeforplots)
        set(gca,'Layer', 'Top')
        
        plot(x_axis_values,occupancy_higher_stored_Unmodified,'w','linewidth',0.5)
       
        h1 = fill([0 x_axis_values, x_axis_values(end)],...        
          [100, occupancy_higher_stored_Unmodified, 100],...
          'w','EdgeColor','none');
            
        xlabel('Time (min)')
        ylabel('Percent Occupancy','fontsize',15)
        ylim([0,100])
        box off
        set(gca,'fontsize',graphfontsize)
    
    subplot(1,4,3)
        plot(x_axis_values, occupancy_higher_stored_phosform1,'w', 'linewidth',0.5)
        hold on
         h1 = fill([0 x_axis_values, x_axis_values(end)],...        
          [1, occupancy_higher_stored_phosform1, 1],...
          'g','EdgeColor','none');
        set(h1,'facealpha',.15)
    
        plot(x_axis_values, occupancy_lower_stored_phosform1,'w','linewidth',0.5)  
        h2 = fill([0, x_axis_values, x_axis_values(end)],...        
      [1, occupancy_lower_stored_phosform1, 1],...
      'w','EdgeColor','none');
        
        plot(x_axis_values, occupancytrend_single,'g-', 'linewidth',1)%,'marker','.','markersize',markersizeforplots)
        set(gca,'Layer', 'Top')   
    
        xlabel('Time (min)')
        ylabel('Percent Occupancy','fontsize',15)
        ylim([0,100])
        box off
        set(gca,'fontsize',graphfontsize)
    
    subplot(1,4,4)
    plot(x_axis_values, occupancy_higher_stored_phosform2,'w', 'linewidth',0.5)
    hold on
     h1 = fill([0 x_axis_values, x_axis_values(end)],...        
          [1, occupancy_higher_stored_phosform2, 1],...
          [1 .5 0],'EdgeColor','none');
        set(h1,'facealpha',.15)
        plot(x_axis_values, occupancy_lower_stored_phosform2,'w','linewidth',0.5)  
         h2 = fill([0, x_axis_values, x_axis_values(end)],...        
            [1, occupancy_lower_stored_phosform2, 1],...
            'w','EdgeColor','none');
   
        plot(x_axis_values, occupancytrend_double,'color',[1 .5 0],'Linestyle','-','linewidth',1)%,'marker','.','markersize',markersizeforplots)
        set(gca,'Layer', 'Top')   

        xlabel('Time (min)')
        ylabel('Percent Occupancy','fontsize',15)
        ylim([0,100])
        box off
        set(gca,'fontsize',graphfontsize)

        xlabel('Time (min)')
        ylabel('Percent Occupancy')
        ylim([0,100])
        box off
        set(gca,'fontsize',graphfontsize)
end

