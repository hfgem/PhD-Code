%Null Model Confirmation Code.
%Contains blocks to test the Poisson and Hexagonal null model density to
%average minimum distance calculations.

%% Poisson Process Null Model

%This section randomly distributes centroids according to a Poisson process
%for different densities and then tests the average minimum distance from
%one point cloud to another to understand the relationship between average
%minimum distance and density.

num_dens = 40;
dens = linspace(0.01,0.99,num_dens); %generate uniform dist vector of densities to test
%size of simulation space
h = 100;
w = 100;
A = h*w;
%generate cloud for each density, calculate weighted average minimum
%distance, and store
avg_min_dist_by_dens = zeros(num_dens,num_dens);
for i = 1:num_dens
    %Find number of points in first cloud to generate
    num_pts1 = floor(dens(i)*A);
    %Generate locations of first density cloud points randomly
    x_pts1 = randi([0,w],[num_pts1,1]);
    y_pts1 = randi([0,h],[num_pts1,1]);
    for j = i:num_dens
        %Find number of points in second cloud to generate from Poisson
        num_pts2 = floor(dens(j)*A);
        %Generate locations of second density cloud points randomly
        x_pts2 = randi([0,w],[num_pts2,1]);
        y_pts2 = randi([0,h],[num_pts2,1]);
        %Calculate the average minimum distance
        x_mat1 = x_pts1.*ones(size(x_pts2))';
        x_mat2 = ones(size(x_pts1)).*x_pts2';
        x_diff = (x_mat1 - x_mat2); %nondim
        y_mat1 = y_pts1.*ones(size(y_pts2))';
        y_mat2 = ones(size(y_pts1)).*y_pts2';
        y_diff = (y_mat1 - y_mat2); %nondim
        all_dist = sqrt(x_diff.^2 + y_diff.^2);
        min_dist1 = min(all_dist,[],2);
        min_dist2 = min(all_dist,[],1);
        w_avg_min_dist = (num_pts1*mean(min_dist1) + num_pts2*mean(min_dist2))/(num_pts1 + num_pts2);
        avg_min_dist_by_dens(i,j) = w_avg_min_dist;
        avg_min_dist_by_dens(j,i) = w_avg_min_dist;
    end
end
clear i num_pts1 x_pts1 y_pts1 j num_pts2 x_pts2 y_pts2 x_mat1 x_mat2 ...
    x_diff y_mat1 y_mat2 y_diff all_dist min_dist1 min_dist2 w_avg_min_dist
%compute factor
factor = avg_min_dist_by_dens.*sqrt(dens);
x = reshape(factor,1,[]);
%histogram
figure;
histogram(x,num_dens)
title('Histogram of Poisson Null Simulation Factors')
xlabel('Factor Value')
ylabel('Number of Density Combinations')
xline(0.5)

%% Hexagonal Distribution Null Model

num_dens = 40;
dens = linspace(0.01,0.99,num_dens);
%size of simulation space
s = 100; %side of simulation space (note: s^2*min(dens) must be >= 1)
%create hexagonal grids
hex_ind = hexagonal_null(s,dens); %struct with x and y ind
%compute avg distances
avg_min_dist_by_dens = zeros(num_dens,num_dens);
for i = 1:num_dens
    x_1 = [hex_ind(i).x_ind];
    y_1 = [hex_ind(i).y_ind];
    n_1 = length(x_1);
    for j = i:num_dens
        x_2 = [hex_ind(j).x_ind];
        y_2 = [hex_ind(j).y_ind];
        n_2 = length(x_2);
        x_mat1 = x_1'*ones(size(x_2));
        y_mat1 = y_1'*ones(size(y_2));
        x_mat2 = ones(size(x_1))'*x_2;
        y_mat2 = ones(size(y_1))'*y_2;
        x_diff = x_mat1 - x_mat2;
        y_diff = y_mat1 - y_mat2;
        all_dist = sqrt(x_diff.^2 + y_diff.^2);
        min_dist1 = min(all_dist,[],2);
        min_dist2 = min(all_dist,[],1);
        w_avg_min_dist = (n_1*mean(min_dist1) + n_2*mean(min_dist2))/(n_1 + n_2);
        avg_min_dist_by_dens(i,j) = w_avg_min_dist;
        avg_min_dist_by_dens(j,i) = w_avg_min_dist;
    end
end
clear i x_1 y_1 n_1 j x_2 y_2 n_2 x_mat1 y_mat1 x_mat2 y_mat2 x_diff ...
    y_diff all_dist min_dist1 min_dist2 w_avg_min_dist
% dist = hexagonal_dist(hex_ind,s,dens);
%compute factor
factor = avg_min_dist_by_dens.*sqrt(dens);
x = reshape(factor,1,[]);
%histogram
figure;
histogram(x,num_dens)
title('Histogram of Hexagonal Null Simulation Factors')
xlabel('Factor Value')
ylabel('Number of Density Combinations')
xline(0.4)
