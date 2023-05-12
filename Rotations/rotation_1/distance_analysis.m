%Analysis of Distance Data from NMJ Macro Outputs
%Any portion labeled "INPUT" requires user input - instructions given

%% Initialization
%INPUT: provide filepath where all of the outputs of the macro exist
filepath = "/Users/hannahgermaine/Documents/PhD/Rotations/Thomas_Fai/PAZ_analysis/AnalysisV9.92fullrandom/";
%INPUT: provide the first portion of the individual image folder names
data_folder_name = "MAX_cleaned";
%automatically get the full list of folders to run through
foldernames = dir(strcat(filepath,strcat(data_folder_name,'*')));

%% Store Centroid Data Per Neuron
centroid_data = struct; %reformat data into a structure
centroid_indices = struct; %store just the centroid indices for each protein
orig_dist = struct; %store the macro distance values (in microns!)

for i = 1:length(foldernames)
    full_folder_name = foldernames(i).name; %name of cell folder
    spltstr = split(full_folder_name,'#'); %grab the cell number
    cell_num = str2num(string(spltstr(2))); %#ok<ST2NM>
    centroid_data(i).cell_num = cell_num; %#ok<*SAGROW> %create new field in structure for cell
    centroid_indices(i).cell_num = cell_num;
    full_folder_path = strcat(filepath,full_folder_name,'/');
    folder_files = dir(strcat(full_folder_path,'*nnd.csv')); %find all centroid loc / distance files
    %Run through all centroid/distance files and store to structure
    for j = 1:length(folder_files)
        dist_filename = folder_files(j).name; %name of file
        %Generate Field Name
        spltstr = split(dist_filename,'--');
        first_prot = string(spltstr(1));
        first_prot = replace(first_prot,'-','_');
        spltstr = split(spltstr(2),'_');
        second_prot = string(spltstr(1));
        second_prot = replace(second_prot,'-','_');
        new_fieldname1 = strcat(first_prot,'_',second_prot);
        new_fieldname2 = strcat(second_prot,'_',first_prot);
        %Grab data
        data = readmatrix(strcat(full_folder_path, dist_filename));
        first_points = data(:,1:2);
        nan1 = isnan(first_points);
        first_points(nan1(:,1),:) = [];
        closest_set_1 = data(~nan1(:,1),1:4);
        dist_1 = data(~nan1(:,1),5);
        second_points = data(:,6:7);
        nan2 = isnan(second_points);
        second_points(nan2(:,1),:) = [];
        closest_set_2 = data(~nan2(:,1),6:9);
        dist_2 = data(~nan2(:,1),10);
        %Store in struct
        centroid_data(i).(first_prot) = first_points;
        centroid_indices(i).(first_prot) = first_points;
        centroid_data(i).(second_prot) = second_points;
        centroid_indices(i).(second_prot) = second_points;
        centroid_data(i).(new_fieldname1) = closest_set_1;
        centroid_data(i).(new_fieldname2) = closest_set_2;
        centroid_data(i).(strcat(first_prot,'_',second_prot,'_mindist')) = dist_1;
        centroid_data(i).(strcat(second_prot,'_',first_prot,'_mindist')) = dist_2;
        orig_dist(i).(strcat(first_prot,'_',second_prot)) = mean(dist_1);
        orig_dist(i).(strcat(second_prot,'_',first_prot)) = mean(dist_2);
    end
end

clear i full_folder_name spltstr cell_num full_folder_path folder_files ...
    j dist_filename first_prot nan1 second_prot nan2 new_fieldname1 new_fieldname2 ...
    data first_points second_points closest_set_1 dist_1 closest_set_2 ...
    dist_2

save(strcat(filepath,'original_dist.mat'),'orig_dist')
save(strcat(filepath,'centroid_data.mat'),'centroid_data')
save(strcat(filepath,'centroid_indices.mat'),'centroid_indices')

%% Calculate Average Min Dist

%INPUT: If data not in path already, uncomment next lines
load(strcat(filepath,'centroid_indices.mat'))

fn = fieldnames(centroid_indices);
prot_num = length(fn)-1;
cell_num = length([centroid_indices.cell_num]);

%Run through each cell and all fields and calculate the average minimum 
%normalized distance

%In the macro, the normalization of distance appears to be based on the
%length of the largest trace. For the sake of this code, because I do not
%have the same data the macro has, I will normalize based on the largest
%centroid-centroid distance in each dimension (imagine there's a centroid
%at each the boundary of a rectangle, and the distance between them is 1).

dist = struct; %to store all data
for i = 1:cell_num
    dist(i).names = {};
    dist(i).avg = zeros(prot_num,prot_num);
    dist(i).std_avg = zeros(prot_num,prot_num);
    dist(i).wavg = zeros(prot_num,prot_num);
    dist(i).counts = zeros(1,prot_num);
    for j = 2:length(fn) %first prot
        for k = j:length(fn) %second prot
            %store name of relationship
            comb_name = [fn{j},'-',fn{k}];
            dist(i).names(j-1,k-1) = {comb_name};
            dist(i).names(k-1,j-1) = {comb_name};
            %grab the centroid indices
            first_cent = [centroid_indices(i).(fn{j})];
            second_cent = [centroid_indices(i).(fn{k})];
            dist(i).counts(1,j-1) = length(first_cent); %store the centroid counts
            dist(i).counts(1,k-1) = length(second_cent); %store the centroid counts
            %grab centroid values
            y1 = first_cent(:,2);
            x1 = first_cent(:,1);
            y2 = second_cent(:,2);
            x2 = second_cent(:,1);
            %grab numbers of centroids
            n1 = length(x1);
            n2 = length(x2);
            dist(i).counts(1,j-1) = n1;
            dist(i).counts(1,k-1) = n2;
            %calculate the y range and x range to use in normalization
            y_norm = max(max(y1),max(y2)) - min(min(y1),min(y2));
            x_norm = max(max(x1),max(x2)) - min(min(x1),min(x2));
            %calculate distances
            x_mat1 = x1*ones([1,n2]);
            x_diff = (x_mat1' - x2)'./x_norm; %normalized dist
            y_mat1 = y1*ones([1,n2]);
            y_diff = (y_mat1'-y2)'./y_norm; %normalized dist
            dist_mat = sqrt(x_diff.^2 + y_diff.^2);
            closest_dist1 = min(dist_mat,[],2); %forward direction
            closest_dist2 = min(dist_mat,[],1); %backward direction
            %store averages in each direction and error values
            dist(i).avg(j-1,k-1) = mean(closest_dist1);
            dist(i).avg(k-1,j-1) = mean(closest_dist2);
            dist(i).std_avg(j-1,k-1) = std(closest_dist1);
            dist(i).std_avg(k-1,j-1) = std(closest_dist2);
            %calculate the weighted average delta and store
            wavg = (n1*mean(closest_dist1) + n2*mean(closest_dist2))/(n1 + n2);
            dist(i).wavg(j-1,k-1) = wavg;
            dist(i).wavg(k-1,j-1) = wavg;
            clear comb_name first_cent second_cent y1 x1 y2 x2 y_norm ...
                x_norm n1 n2 x_mat1 x_diff y_mat1 y_diff dist_mat ...
                closest_dist1 closest_dist2 wavg
        end
    end
    clear j k
end
clear i

save(strcat(filepath,'dist.mat'),'dist')

%% Calculate Factor Vals from W-avg Distances and Compare to Null Models
%INPUT: Uncomment if not in path already
load(strcat(filepath,'dist.mat'))

%Calculate the factor associated with the bootstrapped distance vals
[~,l] = size(dist);
factor = []; %to store all calculated factors across all cells and protein combinations
for cell = 1:l
    dist_mat = dist(cell).wavg; %load the weighted average distances
    counts = dist(cell).counts;
    [num_prot,~] = size(dist_mat);
    cell_factor = [];
    for j = 1:num_prot
        n1 = counts(1,j);
        for k = j+1:num_prot %to avoid the 0 dist of same-protein combinations
            n2 = counts(1,k);
            %the next formula is a reformulation of the density
            %relationship using the number of centroids only since area has
            %been normalized to 1 in the previous section
            dens_weight = (n1^2*sqrt(n2) + n2^2*sqrt(n1))/(n1^2*n2 + n2^2*n1); %density weighted average
            fac = dist_mat(j,k)/dens_weight;
            cell_factor(end+1) = fac;
            ind = ind+1;
        end
    end
    factor(end+1,:) = cell_factor;
end
save(strcat(filepath,'factor_wavg.mat'),'factor')

%Factor Histogram
figure;
x = reshape(factor,1,[]);
histogram(x,ceil(length(x)/10))
title('Calculated Fly NMJ Centroid Factor Distribution')
xlabel('Factor Value')
ylabel('Number of Factors')

%Calculate variation from Poisson factor
p_fac = 0.5;
p_dist = (p_fac - factor)/p_fac;
p_avg = mean(p_dist,'all');
p_std = std(p_dist,[],'all');
figure;
x = 100*reshape(p_dist,1,[]);
histogram(x,ceil(length(x)/10))
title('Variation from Poisson Null Model Factor')
xlabel('Percent Difference')
ylabel('Number of Factors')

%Calculate variation from Hexagonal factor
h_fac = 0.4;
h_dist = (h_fac - factor)/h_fac;
h_avg = mean(h_dist,'all');
h_std = std(h_dist,[],'all');
figure;
x = 100*reshape(h_dist,1,[]);
histogram(x,ceil(length(x)/10))
title('Variation from Hexagonal Null Model Factor')
xlabel('Percent Difference')
ylabel('Number of Factors')

%% Calculate Deltas to Random Points

%INPUT: If data not in path already, uncomment next lines
load(strcat(filepath,'centroid_indices.mat'))

%Basic variables of the centroid data
fn = fieldnames(centroid_indices);
cell_num = length([centroid_indices.cell_num]);

%Run through all cells and proteins and throw random points on a plane to
%see the average minimum distance and calculated factor
%INPUT: number of times to create samples
num_iter = 100;

dist = zeros(cell_num,length(fn)-1); %store avg min dist
err = zeros(cell_num,length(fn)-1); %store error
factor = zeros(cell_num,length(fn)-1); %store calculated factor value

for i = 1:cell_num
    for j = 2:length(fn)
        centroids = centroid_indices(i).(fn{j});
        cent_y = centroids(:,2);
        cent_x = centroids(:,1);
        n1 = length(cent_y);
        y_max = ceil(max(cent_y));
        y_min = floor(min(cent_y));
        x_max = ceil(max(cent_x));
        x_min = floor(min(cent_x));
        y_norm = y_max - y_min;
        x_norm = x_max - x_min;
        %throw iter number of points down to calculate minimum distances
        iter_dist = [];
        for k = 1:num_iter
            rand_x = randi([x_min,x_max]);
            rand_y = randi([y_min,y_max]);
            x_diff = (cent_x - rand_x*ones(size(cent_x)))/x_norm; %normalize so area = 1
            y_diff = (cent_y - rand_y*ones(size(cent_y)))/y_norm; %normalize so area = 1
            dist_pts = sqrt(x_diff.^2 + y_diff.^2);
            iter_dist(end+1) = min(dist_pts);
        end
        dist(i,j-1) = mean(iter_dist); %store avg min dist
        err(i,j-1) = std(iter_dist); %store error of min dist
        %calculate factor
        factor(i,j-1) = dist(i,j-1)*sqrt(n1);
        %since the area is 1, the density is just the number of points, so
        %we have delta = factor/sqrt(num)
    end
end 

clear i j centroids cent_y cent_x y_max y_min x_max x_min boot_dist k ...
    rand_x rand_y x_diff y_diff n dens_weight

if ~isfolder(strcat(filepath,'rand_pt/'))
    mkdir(strcat(filepath,'rand_pt/'));
end

save(strcat(filepath,'rand_pt/','distances.mat'),'dist')
save(strcat(filepath,'rand_pt/','error.mat'),'err')
save(strcat(filepath,'rand_pt/','factor.mat'),'factor')

figure;
x = reshape(factor,1,[]);
histogram(x,ceil(length(x)/6));
title('Distribution of Factors from Centroid to Random Point')
xlabel('Factor Value')
ylabel('Number of Factors')

%% Calculated Deltas to Poisson Null Model Points

%INPUT: If data not in path already, uncomment next lines
load(strcat(filepath,'centroid_indices.mat'))

%Basic variables of the centroid data
fn = fieldnames(centroid_indices);
cell_num = length([centroid_indices.cell_num]);

%Run through all cells and proteins and throw poisson distributed points on
%the centroid plane in different densities to calculate average minimum
%distances and calculate the factor
%INPUT: number of densities to test for the Poisson points
num_dens = 40;
dens = linspace(0.01,0.99,num_dens);

dist = struct; %store avg min dist, error, and factor

for i = 1:cell_num
%     dist(i).('Avg_min_dist_prot') = zeros(length(fn)-1,num_dens); %average min distance averaged over bootstrapped samples
%     dist(i).('Error_prot') = zeros(length(fn)-1,num_dens); %error in the bootstrapped average min dist
%     dist(i).('Factor_prot') = zeros(length(fn)-1,num_dens);
%     dist(i).('Avg_min_dist_poiss') = zeros(num_dens,length(fn)-1); %average min distance averaged over bootstrapped samples
%     dist(i).('Error_poiss') = zeros(num_dens,length(fn)-1); %error in the bootstrapped average min dist
%     dist(i).('Factor_poiss') = zeros(num_dens,length(fn)-1);
    dist(i).('Wavg_min_dist') = zeros(num_dens,length(fn)-1); %weighted average minimum distances
    dist(i).('Factor_wavg') = zeros(num_dens,length(fn)-1); %weighted average factors
    dist(i).('Proteins') = {fn{2:length(fn)}};
    dist(i).('Poisson_Densities') = dens;
    for j = 2:length(fn)
        %grab centroid indices
        centroids = centroid_indices(i).(fn{j});
        cent_y = centroids(:,2);
        cent_x = centroids(:,1);
        num_cent = length(cent_x);
        %calculate norm values
        y_max = ceil(max(cent_y));
        y_min = floor(min(cent_y));
        x_max = ceil(max(cent_x));
        x_min = floor(min(cent_x));
        y_norm = y_max - y_min;
        x_norm = x_max - x_min;
        area = y_norm*x_norm;
        %throw density-determined number of random points down to calculate
        %minimum distances
        for k = 1:num_dens
            num_pts = ceil(dens(k)*area);
            %Create Poisson points
            rand_x = randi([x_min,x_max],1,num_pts);
            rand_y = randi([y_min,y_max],1,num_pts);
            %Calculate normalized distances
            x_mat_1 = cent_x*ones(1,num_pts); %x points in centroids (num_cent x num_pts)
            x_mat_2 = ones(num_cent,1)*rand_x; %x points in Poisson (num_cent x num_pts)
            y_mat_1 = cent_y*ones(1,num_pts); %y points in centroids (num_cent x num_pts)
            y_mat_2 = ones(num_cent,1)*rand_y; %y points in Poisson (num_cent x num_pts)
            x_diff = (x_mat_1 - x_mat_2)/x_norm; %normalize so area = 1
            y_diff = (y_mat_1 - y_mat_2)/y_norm; %normalize so area = 1
            dist_pts = sqrt(x_diff.^2 + y_diff.^2);
            min_1 = min(dist_pts,[],2);
            min_2 = min(dist_pts,[],1);
            %forward dir
            avg_min_1 = mean(min_1);
            dens_weight = (num_cent^2*sqrt(num_pts) + num_pts^2*sqrt(num_cent))/(num_cent^2*num_pts + num_pts^2*num_cent);
%             dist(i).('Avg_min_dist_prot')(j-1,k) = avg_min_1;
%             dist(i).('Error_prot')(j-1,k) = std(min_1);
%             dist(i).('Factor_prot')(j-1,k) = avg_min_1/dens_weight;
            %backward dir
            avg_min_2 = mean(min_2);
%             dist(i).('Avg_min_dist_poiss')(k,j-1) = avg_min_2;
%             dist(i).('Error_poiss')(k,j-1) = std(min_2);
%             dist(i).('Factor_poiss')(k,j-1) = avg_min_2/dens_weight;
            %weighted average
            wavg = (num_cent*avg_min_1 + num_pts*avg_min_2)/(num_cent + num_pts);
            dist(i).('Wavg_min_dist')(j-1,k) = wavg;
            dist(i).('Factor_wavg')(j-1,k) = wavg/dens_weight;
            clear num_pts rand_x rand_y x_mat_1 x_mat_2 y_mat_1 y_mat_2 ...
                x_diff y_diff dist_pts min_1 min_2 avg_min_1 avg_min_2 ...
                dens_weight wavg
        end
        clear centroids cent_y cent_x num_cent y_max y_min x_max x_min ...
            y_norm x_norm area
    end
end

clear i j

if ~isfolder(strcat(filepath,'poisson_null/'))
    mkdir(strcat(filepath,'poisson_null/'));
end

save(strcat(filepath,'poisson_null/','dist_regular.mat'),'dist')

% avg_avg_min_dist_prot = mean([dist.Avg_min_dist_prot],'all');
% avg_avg_min_dist_poiss = mean([dist.Avg_min_dist_poiss],'all');
avg_wavg_min_dist = mean([dist.Wavg_min_dist],'all');

figure;
% a = reshape([dist.Avg_min_dist_prot],1,[]);
% b = reshape([dist.Avg_min_dist_poiss],1,[]);
c = reshape([dist.Wavg_min_dist],1,[]);
% agrp = ones(size(a));
% bgrp = 2*ones(size(b));
cgrp = 3*ones(size(c));
%boxplot([a,b,c],[agrp,bgrp,cgrp],'Labels',{'Centroid to Poisson','Poisson to Centroid','Weighted Average'})
boxplot(c,cgrp)
%title('Poisson Deltas')
title('Centroid to Poisson Deltas')

figure;
% histogram(a,length(a)/10);
% hold on
% histogram(b,length(b)/10);
histogram(c,ceil(sqrt(length(c))));
title('Distribution of Factors from Centroid to Spatial Poisson Distribution')
xlabel('Factor Value')
ylabel('Number of Factors')
%legend('Centroid to Poisson','Poisson to Centroid','Weighted Average')

%% Calculated Deltas to Poisson Null Model Points with Bootstrapping
%NEEDS UPDATING WITH WAVG CALC

%INPUT: If data not in path already, uncomment next lines
load(strcat(filepath,'centroid_indices.mat'))

%Basic variables of the centroid data
fn = fieldnames(centroid_indices);
cell_num = length([centroid_indices.cell_num]);

%Run through all cells and proteins and throw poisson distributed points on
%the centroid plane in different densities to calculate average minimum
%distances and calculate the factor
%INPUT: bootstrap iteration number
num_boot = 50;
%INPUT: fraction of dataset to sample each iteration
samp_frac = 0.75;
%INPUT: number of densities to test for the Poisson points
num_dens = 40;
dens = linspace(0.01,0.99,num_dens);

dist = struct; %store avg min dist, error, and factor

for i = 1:cell_num
%     dist(i).('Avg_min_dist_prot') = zeros(length(fn)-1,num_dens); %average min distance from poisson to protein
%     dist(i).('Avg_min_dist_poiss') = zeros(num_dens,length(fn)-1); %average min distance from protein to poisson
%     dist(i).('Error_prot') = zeros(length(fn)-1,num_dens); %error in the bootstrapped average min dist
%     dist(i).('Error_poiss') = zeros(num_dens,length(fn)-1); %error in the bootstrapped average min dist
%     dist(i).('Factor_prot') = zeros(length(fn)-1,num_dens);
%     dist(i).('Factor_poiss') = zeros(num_dens,length(fn)-1);
    dist(i).('Wavg_min_dist') = zeros(num_dens,length(fn)-1); %weighted average minimum distances
    dist(i).('Factor_wavg') = zeros(num_dens,length(fn)-1); %weighted average factors
    dist(i).('Proteins') = {fn{2:length(fn)}};
    dist(i).('Poisson_Densities') = dens;
    dist(i).('Poisson_Pt_Counts') = zeros(length(fn)-1,num_dens);
    for j = 2:length(fn)
        centroids = centroid_indices(i).(fn{j});
        cent_y = centroids(:,2);
        cent_x = centroids(:,1);
        num_cent = length(cent_x);
        y_max = ceil(max(cent_y));
        y_min = floor(min(cent_y));
        x_max = ceil(max(cent_x));
        x_min = floor(min(cent_x));
        y_norm = y_max - y_min;
        x_norm = x_max - x_min;
        area = y_norm*x_norm;
        %throw density-proper number of points down to calculate minimum distances
        for k = 1:num_dens
            num_pts = ceil(dens(k)*area);
            dist(i).('Poisson_Pt_Counts')(j-1,k) = num_pts;
            boot_avg_min_1 = [];
            boot_avg_min_2 = [];
            boot_wavg = [];
            %Create Poisson points
            rand_x = randi([x_min,x_max],1,num_pts);
            rand_y = randi([y_min,y_max],1,num_pts);
            %Bootstrap avg min dist calculation from Poisson points
            for l = 1:num_boot
                boot_ind = randi([1,floor(samp_frac*num_boot)]);
                rand_x_boot = rand_x(boot_ind);
                rand_y_boot = rand_y(boot_ind);
                %Calculate normalized distances
                x_mat_1 = cent_x*ones(1,num_pts); %x points in centroids (num_cent x num_pts)
                x_mat_2 = ones(num_cent,1)*rand_x; %x points in Poisson (num_cent x num_pts)
                y_mat_1 = cent_y*ones(1,num_pts); %y points in centroids (num_cent x num_pts)
                y_mat_2 = ones(num_cent,1)*rand_y; %y points in Poisson (num_cent x num_pts)
                x_diff = (x_mat_1 - x_mat_2)/x_norm; %normalize so area = 1
                y_diff = (y_mat_1 - y_mat_2)/y_norm; %normalize so area = 1
                dist_pts = sqrt(x_diff.^2 + y_diff.^2);
                min_1 = min(dist_pts,[],2);
                boot_avg_min_1(end+1) = mean(min_1);
                min_2 = min(dist_pts,[],1);
                boot_avg_min_2(end+1) = mean(min_2);
                wavg = (num_cent*min_1 + num_pts*min_2)/(num_cent + num_pts);
                boot_wavg = mean(wavg);
            end
            %Store bootstrapped averages and vals
            %forward dir
%             avg_min_1 = mean(boot_avg_min_1);
            dens_weight = (num_cent^2*sqrt(num_pts) + num_pts^2*sqrt(num_cent))/(num_cent^2*num_pts + num_pts^2*num_cent);
%             dist(i).('Avg_min_dist_prot')(j-1,k) = avg_min_1;
%             dist(i).('Error_prot')(j-1,k) = std(min(dist_pts,[],2));
%             dist(i).('Factor_prot')(j-1,k) = avg_min_1/dens_weight;
            %backward dir
%             avg_min_2 = mean(boot_avg_min_2);
%             dist(i).('Avg_min_dist_poiss')(k,j-1) = avg_min_2;
%             dist(i).('Error_poiss')(k,j-1) = std(min(dist_pts,[],1));
%             dist(i).('Factor_poiss')(k,j-1) = avg_min_2/dens_weight;
            %weighted average
            wavg = mean(boot_wavg);
            dist(i).('Wavg_min_dist')(j-1,k) = wavg;
            dist(i).('Factor_wavg')(j-1,k) = wavg/dens_weight;
        end
    end
end 

clear i j centroids cent_y cent_x y_max y_min x_max x_min boot_dist k ...
    rand_x rand_y x_diff y_diff n dens_weight avg_min_1 avg_min_2 x_mat_1 ...
    x_mat_2 y_mat_1 y_mat_2 area num_pts num_dens num_cent min_1 min_2 ...
    boot_avg_min_1 boot_avg_min_2 l x_norm y_norm dist_pts

if ~isfolder(strcat(filepath,'poisson_null/'))
    mkdir(strcat(filepath,'poisson_null/'));
end

save(strcat(filepath,'poisson_null/','dist_boot.mat'),'dist')

figure;
% a = reshape([dist.Avg_min_dist_prot],1,[]);
% b = reshape([dist.Avg_min_dist_poiss],1,[]);
c = reshape([dist.Wavg_min_dist],1,[]);
% agrp = ones(size(a));
% bgrp = 2*ones(size(b));
cgrp = 3*ones(size(c));
%boxplot([a,b,c],[agrp,bgrp,cgrp],'Labels',{'Centroid to Poisson','Poisson to Centroid','Weighted Average'})
boxplot(c,cgrp)
%title('Poisson Deltas')
title('Centroid to Poisson Deltas with Bootstrapping')

figure;
% histogram(a,length(a)/10);
% hold on
% histogram(b,length(b)/10);
histogram(c,ceil(sqrt(length(c))));
title('Distribution of Factors from Centroid to Spatial Poisson Distribution with Bootstrapping')
xlabel('Factor Value')
ylabel('Number of Factors')
%legend('Centroid to Poisson','Poisson to Centroid','Weighted Average')

%% Calculate Deltas to Hexagonal Null Model Points

%INPUT: If data not in path already, uncomment next lines
load(strcat(filepath,'centroid_indices.mat'))

%Basic variables of the centroid data
fn = fieldnames(centroid_indices);
cell_num = length([centroid_indices.cell_num]);

%Run through all cells and proteins and throw poisson distributed points on
%the centroid plane in different densities to calculate average minimum
%distances and calculate the factor
%INPUT: number of densities to test for the Hexagonal lattice
num_dens = 40;
dens = linspace(0.01,0.99,num_dens);

dist = struct; %store avg min dist, error, and factor

for i = 1:cell_num
%     dist(i).('Avg_min_dist_prot') = zeros(length(fn)-1,num_dens); %average min distance averaged over bootstrapped samples
%     dist(i).('Error_prot') = zeros(length(fn)-1,num_dens); %error in the bootstrapped average min dist
%     dist(i).('Factor_prot') = zeros(length(fn)-1,num_dens);
%     dist(i).('Avg_min_dist_hex') = zeros(num_dens,length(fn)-1); %average min distance averaged over bootstrapped samples
%     dist(i).('Error_hex') = zeros(num_dens,length(fn)-1); %error in the bootstrapped average min dist
%     dist(i).('Factor_hex') = zeros(num_dens,length(fn)-1);
    dist(i).('Wavg_min_dist') = zeros(length(fn)-1,num_dens); %weighted average minimum distances
    dist(i).('Factor_wavg') = zeros(length(fn)-1,num_dens); %weighted average factor values
    dist(i).('Proteins') = {fn{2:length(fn)}};
    dist(i).('Hexagonal_Densities') = dens;
    for j = 2:length(fn)
        %grab centroid indices
        centroids = centroid_indices(i).(fn{j});
        cent_y = centroids(:,2);
        cent_x = centroids(:,1);
        num_cent = length(cent_x);
        %calculate norm values
        y_max = ceil(max(cent_y));
        y_min = floor(min(cent_y));
        x_max = ceil(max(cent_x));
        x_min = floor(min(cent_x));
        y_norm = y_max - y_min;
        x_norm = x_max - x_min;
        %INPUT: MAKE SURE THE hexagonal_null_2 FUNCTION IS IN THE DIRECTORY
        hex_ind = hexagonal_null_2(y_norm,x_norm,dens); %indices of Hex lattices
        %throw boot number of points down to calculate minimum distances
        for k = 1:num_dens
            %Grab Hexagonal lattice points
            hex_x = [hex_ind(k).x_ind] + x_min;
            hex_y = [hex_ind(k).y_ind] + y_min;
            num_pts = ceil(length(hex_x));
            %Calculate normalized distances
            x_mat_1 = cent_x*ones(1,num_pts); %x points in centroids (num_cent x num_pts)
            x_mat_2 = ones(num_cent,1)*hex_x; %x points in Poisson (num_cent x num_pts)
            y_mat_1 = cent_y*ones(1,num_pts); %y points in centroids (num_cent x num_pts)
            y_mat_2 = ones(num_cent,1)*hex_y; %y points in Poisson (num_cent x num_pts)
            x_diff = (x_mat_1 - x_mat_2)/x_norm; %normalize so area = 1
            y_diff = (y_mat_1 - y_mat_2)/y_norm; %normalize so area = 1
            dist_pts = sqrt(x_diff.^2 + y_diff.^2);
            min_1 = min(dist_pts,[],2);
            min_2 = min(dist_pts,[],1);
            %forward dir
            avg_min_1 = mean(min_1);
            dens_weight = (num_cent^2*sqrt(num_pts) + num_pts^2*sqrt(num_cent))/(num_cent^2*num_pts + num_pts^2*num_cent);
%             dist(i).('Avg_min_dist_prot')(j-1,k) = avg_min_1;
%             dist(i).('Error_prot')(j-1,k) = std(min_1);
%             dist(i).('Factor_prot')(j-1,k) = avg_min_1/dens_weight;
            %backward dir
            avg_min_2 = mean(min_2);
%             dens_weight = (num_pts^2*sqrt(num_cent) + num_cent^2*sqrt(num_pts))/(num_pts^2*num_cent + num_cent^2*num_pts);
%             dist(i).('Avg_min_dist_hex')(k,j-1) = avg_min_2;
%             dist(i).('Error_hex')(k,j-1) = std(min_2);
%             dist(i).('Factor_hex')(k,j-1) = avg_min_2/dens_weight;
            %weighted average
            wavg = (num_cent*avg_min_1 + num_pts*avg_min_2)/(num_cent + num_pts);
            dist(i).('Wavg_min_dist')(j-1,k) = wavg;
            dist(i).('Factor_wavg')(j-1,k) = wavg/dens_weight;
        end
    end
end 

clear i j centroids cent_y cent_x y_max y_min x_max x_min boot_dist k ...
    hex_x hex_y x_diff y_diff n dens_weight avg_min_1 avg_min_2 x_mat_1 ...
    x_mat_2 y_mat_1 y_mat_2 area num_pts num_dens num_cent min_1 min_2 ...
    x_norm y_norm dist_pts

if ~isfolder(strcat(filepath,'hex_null/'))
    mkdir(strcat(filepath,'hex_null/'));
end

save(strcat(filepath,'hex_null/','dist_regular.mat'),'dist')

figure;
% a = reshape([dist.Avg_min_dist_prot],1,[]);
% b = reshape([dist.Avg_min_dist_hex],1,[]);
c = reshape([dist.Wavg_min_dist],1,[]);
% agrp = ones(size(a));
% bgrp = 2*ones(size(b));
cgrp = 2*ones(size(c));
%boxplot([a,b],[agrp,bgrp],'Labels',{'Centroid to Hex','Hex to Centroid'})
boxplot(c,cgrp)
title('Centroid to Hexagonal Deltas')

figure;
% histogram(a,ceil(sqrt(length(a))));
% hold on
% histogram(b,ceil(sqrt(length(b))));
histogram(c,ceil(sqrt(length(c))));
title('Distribution of Factors from Centroid to Hexagonal Distribution')
xlabel('Factor Value')
ylabel('Number of Factors')
%legend('Centroid to Hexagonal','Hexagonal to Centroid','Weighted Average')