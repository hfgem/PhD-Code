%Create a base Hexagonal null model and return indices

function hex_ind = hexagonal_null(s,dens)
    A = s^2; %area of space
    hex_ind = struct;
    %For each density create a hexagonal grid
    for i = 1:length(dens)
        a = sqrt(4/(3*sqrt(3)*dens(i))); %edge-length of hexagon
        h = sqrt(3)*a/2; %height of half-hexagon (row height)
        rows = [0:h:s];
        num_rows = length(rows); %how many rows will fit in grid
        %set hexagon indc
        odd_row = [0:a:s];
        odd_clean = [2:3:length(odd_row)];
        odd_row(:,odd_clean) = [];
        even_row = [a/2:a:s-(a/2)];
        even_clean = [3:3:length(even_row)];
        even_row(:,even_clean) = [];
        clear odd_clean even_clean
        %store plotting indices
        x_ind = [];
        y_ind = [];
        for j = 1:num_rows
            y_val = rows(j);
            if mod(j,2) == 1 %oddrow
                y_vec = y_val*ones(size(odd_row));
                x_ind = [x_ind, odd_row];
                y_ind = [y_ind, y_vec];
            elseif mod(j,2) == 0 %evenrow
                y_vec = y_val*ones(size(even_row));
                x_ind = [x_ind, even_row];
                y_ind = [y_ind, y_vec];
            end    
        end
        %Sanity Check
%         figure;
%         set(gca,'Units','normalized','Position',[0 0 1 1]);
%         set(gcf,'Units','pixels','Position',[200 200 2*s 2*s]);
%         scatter(x_ind,y_ind,'.')
        %Store
        hex_ind(i).x_ind = x_ind;
        hex_ind(i).y_ind = y_ind;
    end
end