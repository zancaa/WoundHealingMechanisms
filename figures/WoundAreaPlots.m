%% LOAD DATA AND PLOT
% Need to be in the directory where the results are saved, e.g. the
% WoundMechanisms directory.
% Six simulations types to deal with - note that proliferation is a bit
% different: multiple runs. So is global purse string: wound area is area 
% of last cell
for i = 1:6
    if i == 1
        % Compression
        folder_name = 'Circle/Compression/Cell_target_area_';
        parameter_values = 0.9:0.025:1;
        title_string = 'Compression';
        legend_string = '$A_{\mathrm{target}}=$';
        constant_param = sqrt(3)/2;
    elseif i == 2
        % Differential adhesion (local purse string)
        folder_name = 'Circle/EdgeContractility/CellBoundaryAdhesion_';
        parameter_values = 4:4:20;
        title_string = 'Local purse string';
        legend_string = '$\gamma_{\mathrm{cell,boundary}}=$';
        constant_param = 0;
    elseif i == 3
        % Wound cell (global purse string)
        folder_name = 'Circle/PurseString/WoundCellMembraneEnergy_';
        parameter_values = 0.4:0.4:2;
        title_string = 'Global purse string';
        legend_string = '$\beta_{W}=$';
        constant_param = 0;
    elseif i == 4
        % Normal force (local crawling)
        folder_name = 'Circle/NormalForce/NormalForceStrength_';
        parameter_values = 2:2:10;
        title_string = 'Local crawling';
        legend_string = '$f_{\mathrm{LC}}=$';
        constant_param = 0;
    elseif i == 5
        % Wound centre force (global crawling)
        folder_name = 'Circle/WoundCentreForce/WoundCentreForceStrength_';
        parameter_values = 1:5;
        title_string = 'Global crawling';
        legend_string = '$f_{\mathrm{GC}}=$';
        constant_param = 0;
    elseif i == 6
        % Proliferation
        folder_name = 'Circle/Proliferation/Division_probability_';
        parameter_values = 0.006:0.006:0.3;
        title_string = 'Proliferation';
        legend_string = '$p_{\mathrm{div}}=$';
        constant_param = 0;
    else
        % Do nothing
        return
    end
    
    figure
    hold on
    xlabel('Time (hrs)','Interpreter','latex')
    ylabel('Wound area (CD$^{2}$)','Interpreter','latex')
    title(title_string)
    p0 = plot(0:0.1:20,ones(201,1),'k-','LineWidth',2);
    axis([0 10 0 1.1])
    legend_values = cell(6,1);
    legend_values{1} = strcat(legend_string,num2str(constant_param));
    if i ~= 6  && i ~= 3     
        for j = 1:5
            data_filename = strcat(folder_name,num2str(parameter_values(j)),'/voidArea.dat');
            all_data = importdata(data_filename);
            time = all_data.data(:,1);
            area = all_data.data(:,3);
            area_scaled = area/area(1);
            plot(time,area_scaled,'LineWidth',2)
            legend_values{j+1} = strcat(legend_string,num2str(parameter_values(j)));
        end
        plot(3*ones(length(0:0.1:1.3),1),0:0.1:1.3,'LineStyle','--','LineWidth',2)
        legend(legend_values,'Interpreter','latex','EdgeColor','None','NumColumns',2)
    elseif i == 3
        time = 0:0.1:20;
        area = zeros(201,1);
        for j = 1:5
            data_filename = strcat(folder_name,num2str(parameter_values(j)),...
                '/results_from_time_0/cellareas.dat');
            all_data = importdata(data_filename);
            for k = 1:length(all_data)
                tmp = str2num(all_data{k});
                area(k) = tmp(end);
            end
            area_scaled = area/area(1);
            plot(time,area_scaled,'LineWidth',2)
            legend_values{j+1} = strcat(legend_string,num2str(parameter_values(j)));
        end
        plot(3*ones(length(0:0.1:1.3),1),0:0.1:1.3,'LineStyle','--','LineWidth',2)
        legend(legend_values,'Interpreter','latex','EdgeColor','None','NumColumns',2)
    elseif i == 6
        p = cell(5,1);
        for j = 1:5
            all_runs = zeros(201,10);
            all_runs_scaled = zeros(201,10);
            for k = 1:10
                data_filename = strcat(folder_name,num2str(parameter_values(j)),...
                    '/Run_',num2str(k-1),'/voidArea.dat');
                all_data = importdata(data_filename);
                if length(all_data.data(:,3)) == 201
                    all_runs(:,k) = all_data.data(:,3);
                    all_runs_scaled(:,k) = all_runs(:,k)/all_runs(1,k);
                else
                    tmp = all_data.data(:,3);
                    tmp = [tmp; NaN*ones(201-length(all_data.data(:,3)),1)];
                    all_runs(:,k) = tmp;
                    all_runs_scaled(:,k) = all_runs(:,k)/all_runs(1,k);
                end
            end
            fill([0:0.1:20 fliplr(0:0.1:20)],[max(all_runs_scaled,[],2)' fliplr(min(all_runs_scaled,[],2)')],...
                'k','FaceAlpha',0.3,'EdgeColor','None')
            p{j} = plot(0:0.1:20,median(all_runs_scaled,2),'LineWidth',2);
            legend_values{j+1} = strcat(legend_string,num2str(parameter_values(j)));
        end
        plot(3*ones(length(0:0.1:1.3),1),0:0.1:1.3,'LineStyle','--','LineWidth',2)
        legend([p0 p{1} p{2} p{3} p{4} p{5}],legend_values,'Interpreter','latex','EdgeColor','None','NumColumns',2)
    else
        % Should never be able to reach this
        return
    end
    
end

