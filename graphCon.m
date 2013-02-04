%clear
%clc
%load mytest
%comp = mytest;

function [res tmp_path] = graphCon(comp)
% GRAPH_CON Extract the major skeleton from a given connected component
%
% graphCon(comp) Extract the major skeleton for a given connected component
% comp
%
% res = graphCom(comp) returns a binary matrix who is of the same size as
% comp.
%
% Input:
%       comp :  A binary (0's or 1's) matrix which contains only one
%       connected component. comp should be only one-pixel width.
% Output:
%       res : A binary (0's or 1's) matrix which contains the main skeleton
%       corresponding to comp.
% Example:
%   load mytest;
%   res = graphCon(comp);

% 2009-12-08        Yuan Yuan <yy46@njit.edu>
%tic
[rows,cols] = find(comp);
num = length(rows);
mat = zeros(num);
for i = 1:num
    for j = 1:num
        if i ~= j
            if ( abs(rows(i)-rows(j)) <=1 && abs(cols(i)-cols(j)) <=1 )
%                 if (rows(i) ~= rows(j) && cols(i) ~= cols(j))
%                     mat(i,j) = 1.4;
%                 else
                    mat(i,j) = 1;
%                end
            end
        end
    end
end


%endpoints = find(sum(mat)==1);
%num_end = length(endpoints);
%if num_end >= 2
    
    
    mat2 = sparse(logical(mat));
    mat3 = mst_prim(mat2);
    endpoints = find(sum(mat3)==1);
    num_end = length(endpoints);
    selected_d = [];
    %selected_dt = [];
    selected_pred = [];
    maxd = 0;
    start_point = 0;
    for i = 1:num_end
        %[d pred]=dijkstra(mat2,endpoints(i));
        [d dt pred] = bfs(mat2,endpoints(i));
        tmp = max(max(d));
        if (tmp > maxd)
            maxd = tmp;
            selected_d = d;
            %    selected_dt = dt;
            selected_pred = pred;
            start_point = endpoints(i);
        end
    end
    rst = find( selected_d == max(selected_d),1,'first');
    
    % Print the path
    %fprintf('Path:\n');
    path =[]; u = rst;
    while (u ~= start_point) path=[u path]; u=selected_pred(u); end
    path = [start_point path];
    %fprintf('%s',labels{lax});
    %fprintf('%s',start_point);
    %for i=path; fprintf(' --> %s', labels{i}); end, fprintf('\n');
    %path;
    tmp_path =[];
    res = zeros(size(comp));
    for i = 1:length(path)
        res(rows(path(i)),cols(path(i))) = 1;
        tmp_path = [tmp_path; rows(path(i)),cols(path(i))];
    end
    %pathcell = cell(1,1);
    %path_cell{1} = tmp_path;
%else
%    res = comp;
%end
%toc
