%Citation splitter for overleaf
%Takes citations copy and pasted from Microsoft Word bibliography to Excel
%and saved as citations.csv. NOTE: first value in .csv must be something
%other than a citation. MATLAB readtable otherwise will set the first
%citation as the column name.

%% Import citations

file_path = uigetdir('/Users/hannahgermaine/Documents/PhD/','Select Citation .csv File Folder'); 
citations = readtable(strcat(file_path,'/citations.csv'),'RowNamesColumn',0,'VariableNamesLine',0);

%% Split into LaTeX bib format

%EXAMPLE Format:
% @article{azizi,
%   author =       "Azizi, Amir H. and Wiskott, Laurenz and Cheng, Sen",
%   title =        "A computational model for preplay in the hippocampus",
%   journal =      "Frontiers in Computational Neuroscience",
%   pages =        "1-15",
%   year =         "2013"
% }

bibtext = {};

for i = 1:height(citations)
    orig_cit = citations(i,1).('Var1'){1}; %characters only
    split1 = split(orig_cit,'. (');
    authors = strcat(split1{1},'.');
    author_split = split(authors,'& ');
    if length(author_split) > 1
        authors = [author_split{1} author_split{2}];
    end
    author_split = split(authors,',');
    authors_reworked = '';
    for a = 1:ceil(length(author_split)/2)
        if a < ceil(length(author_split)/2)
            authors_reworked = [authors_reworked author_split{(a-1)*2+1} ',' author_split{a*2} ' and'];
        else
            try
                authors_reworked = [authors_reworked author_split{(a-1)*2+1} ',' author_split{a*2}];
            catch
                authors_reworked = [authors_reworked author_split{(a-1)*2+1}];
            end
        end
    end    
    first_auth = split(authors,',');
    first_auth = first_auth{1};
    split2 = split([split1{2:end}],'). ');
    year = split2{1};
    if length(split2) > 2
        split2{2} = [split2{2},'). '];
    end
    split3 = split([split2{2:end}],'. ');
    if length(split3) == 1
        split3 = split([split2{2:end}],'? ');
        if length(split3) == 1
            split3 = split([split2{2:end}],'! ');
        end
    end    
    if length(split3) > 1
        title = split3{1};
        split4 = split([split3{2:end}],', ');
        if length(split4) == 1
            split4 = split([split3{2:end}],'https://');
            if length(split4) > 1
                url = ['https://',split4{2}];
            end    
        else
            journal = split4{1};
            split5 = split(split4{2},'.');
            pages = split5{1};
        end
    else
        split4 = split(split3,'https://');
        if length(split4) > 1
            title = authors;
            clear authors
            url = ['https://',split4{2}];
        end    
    end
    ae = exist('authors');
    te = exist('title');
    ue = exist('url');
    je = exist('journal');
    pe = exist('pages');
    ye = exist('year');
    if ae == 1
        new_cit = ['@article{' first_auth ',' newline 'author = "' authors '",'];
        if te == 1
             new_cit = [new_cit newline 'title = "' title '"'];
        end
    elseif ue == 1
        if te == 1
            new_cit = ['@online{' title(1:5)];
        else
            spliturl = split(url,'.');
            new_cit = ['@online{' spliturl{2}];
        end
    end
    if je  == 1
        new_cit = [new_cit ',' newline 'journal = "' journal '"'];
    end
    if pe == 1
        new_cit = [new_cit ',' newline 'pages = "' pages '"'];
    end
    if ye == 1
        new_cit = [new_cit ',' newline 'year = "' year '"'];
    end
    if ue == 1
        new_cit = [new_cit ',' newline 'url = "' url '"'];
    end
    bibtext{i} = [new_cit newline '}' newline];
    
    clear orig_cit split1 authors author_split authors_reworked a ...
        first_auth split2 year split3 title split4 url journal pages ...
        je pe ye ue new_cit
end    

writecell(bibtext,strcat(file_path,'/citations.txt'),'QuoteStrings',false,'Delimiter',' ')