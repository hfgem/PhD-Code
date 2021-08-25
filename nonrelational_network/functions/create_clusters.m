function [cluster_mat, conns] = create_clusters(n, clusters, ...
    cluster_n, cluster_prob)
    %_________
    %ABOUT: This function generates the network clusters and connections
    %based on the number of neurons, number of clusters, number of neurons
    %per cluster, and the probability of connections within a cluster.
    %
    %INPUTS:
    %   n = number of neurons in the network
    %   clusters = number of clusters of neurons
    %   cluster_n = number of neurons in a cluster
    %   cluster_prob = intra-cluster connection probability
    %
    %OUTPUTS:
    %   cluster_mat = A binary [clusters x n] matrix of which neurons are in
    %                 which cluster
    %   conns = An [n x n] matrix of which neurons are connected to each
    %                 other, with values greater than 1 implying stronger
    %                 connectivity
    %_________

    %Create clusters by randomly selecting cluster_n neurons for each
    cluster_mat = zeros(clusters,n);
    for i = 1:clusters %set clusters
        ord = randperm(n,cluster_n); %random connectivity
        cluster_mat(i,ord) = 1; %note which neurons are in a cluster
    end
    clear ord i

    %Add back in those neurons that are not placed in a cluster, by
    %removing a place from another neuron with a high presence - this
    %section can be removed if you'd like some neurons to be unconnected
    ind_non = find(sum(cluster_mat) == 0);
    ind_high = find(sum(cluster_mat) > 2);
    for i = ind_non
        clust_place = randi(clusters);
        ind_inclust = find(cluster_mat(clust_place,:));
        val_steal = datasample(intersect(ind_inclust,ind_high),1);
        cluster_mat(clust_place,i) = 1;
        cluster_mat(clust_place,val_steal) = 0;
        ind_high = find(sum(cluster_mat) > 2);
    end
    
    %Find the matrix of total connectivity - it will have integer values of
    %1 or more, representing the strength of connection between neurons
    conns = zeros(n);
    for i = 1:clusters
        ord = find(cluster_mat(i,:));
        ord_len = length(ord);
        conns(ord,ord) = conns(ord,ord) + (rand(ord_len,ord_len) <= cluster_prob);
    end
    
    %Remove self-connectivity
    for i = 1:n
        conns(i,i) = 0;
    end    
end