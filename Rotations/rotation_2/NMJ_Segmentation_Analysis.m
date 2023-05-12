%This program analyzes data output from the PAZ Analysis ImageJ Macro. The
%data must be in a .csv format and the first column should be changed from
%file name to the type (control or experimental) to simplify the analysis.
%Change the column header to 'Type' as well.

%% Data Import and Basic Variables

%User selects the proper .csv file and imports to the table 'data'
[filename, pathname] = uigetfile('*.csv');
path_slash = find(pathname == '/');
savepath = pathname(1:path_slash(end-1));
opts = detectImportOptions(fullfile(pathname,filename));
opts.PreserveVariableNames = 1;
data = readtable(fullfile(pathname,filename),opts);
clear opts path_slash

%Basic Data Info
[rows, columns] = size(data);
headers = data.Properties.VariableNames;

%Collect row indices of control and experimental
type = data(:,1);
control = reshape(find(cellfun(@(x) isequal(x,'Control'),data.Type)),1,[]);
experimental = reshape(find(cellfun(@(x) isequal(x,'Experimental'),data.Type)),1,[]);

%% Gaussian and Same Distribution Analysis

%Specify Data of Interest
%Save the names of the measurements as they appear in the data table before
%the _ (protein name)
measurements = {'NMJ-Mean', 'NMJ-CoV', 'NMJ-PCC', 'Nwk-Cell-Mean', 'Dap160-Cell-Mean', 'BRP-Cell-Mean'};

%Create a structure to store results
results = struct;

for m = measurements
    %Create location in structure to store result
    results.(strrep(m{1},'-','_')) = [];
    
    %Plot box plots for each channel
    col_ind = reshape(find(cellfun(@(x) contains(x,m{1}),headers)),1,[]);
    proteins = {};
    for i = col_ind
        f = figure;
        hold on
        colname = headers{i};
        colnamesplit = split(colname,'_');
        proteins{end+1} = colnamesplit{2};
        title([m{1},' for ',colnamesplit{2}])
        boxplot(data.(colname),data.('Type'))
        set(gca,'fontsize',14)
        saveas(f,strcat(savepath,m{1},'_',colnamesplit{2},'.fig'))
        saveas(f,strcat(savepath,m{1},'_',colnamesplit{2},'.jpg'))
    end
    clear i colname colnamesplit f
    close all

    %Test for Gaussian Distribution
    gaussian = 1; %1 if from normal, 0 if not

    for i = col_ind %if at least one set of data is non-Gaussian, we will treat them all as non-Gaussian
        h_exp = adtest(data{experimental,i});
        h_cont = adtest(data{control,i});
        if h_exp == 1
            gaussian = 0;
            break
        end
        if h_cont == 1
            gaussian = 0;
            break
        end
    end
    results.(strrep(m{1},'-','_')).gaussian_anderson_darling = gaussian;
    clear i h_exp h_cont

    %Run analysis based on distribution type
    if gaussian == 1
        for i = 1:length(col_ind) %for each protein
            [h,p] = ttest2(data{experimental,col_ind(i)},data{control,col_ind(i)}); %0.05 sig level
            results.(strrep(m{1},'-','_')).(strrep(proteins{i},'-','_')) = [h,p];
        end
    else
        for i = 1:length(col_ind) %for each protein
            [p,h] = ranksum(data{experimental,col_ind(i)},data{control,col_ind(i)}); %0.05 sig level
            results.(strrep(m{1},'-','_')).(strrep(proteins{i},'-','_')) = [h,p];
        end
    end
    clear i h p
    
end
clear m

save(strcat(savepath,'results.mat'),'results')

%% Creating Combined Histogram Visualization Across Proteins

measurements = {'NMJ-Mean', 'NMJ-CoV', 'NMJ-PCC', 'Nwk-Cell-Mean', 'Dap160-Cell-Mean', 'BRP-Cell-Mean'};

%Plot box plots for each channel
for m = measurements
    %Collect Data
    col_ind = reshape(find(cellfun(@(x) contains(x,m{1}),headers)),1,[]);
    hist_data = [];
    hist_groups = {};
    for i = col_ind
        colname = headers{i};
        colnamesplit = split(colname,'_');
        proteins{end+1} = colnamesplit{2};
        for j = 1:rows
            hist_data(end+1) = data.(colname)(j);
            hist_groups{end+1} = strcat(colnamesplit{2},'-',data.('Type'){j});
        end
    end
    %Plot Combined Histogram
    f = figure;    
    title(m{1})
    set(gcf,'Position',[0,0,1010,520])
    hold on
    boxplot(hist_data, string(hist_groups))
    set(gca,'fontsize',14)
    saveas(f,strcat(savepath,m{1},'-combined.fig'))
    saveas(f,strcat(savepath,m{1},'-combined.jpg'))
end
