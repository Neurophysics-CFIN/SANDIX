function pars_best = landscape(y,S_fun,grid_points,batch_size,labels)
% grid_points is a cell array with an element for each parameter. Each
% element is a vector of points.

dummy_labels = "$x_" + num2str((1:numel(grid_points))') + "$";
if ~exist('labels','var') || isempty(labels)
    labels = dummy_labels;
end
if numel(labels)<numel(grid_points)
    labels(end+1:numel(grid_points)) = dummy_labels(numel(labels)+1:end);
end

N_points = cellfun(@(x)length(x), grid_points);
if ~exist('batch_size','var') || isempty(batch_size) 
    batch_size = prod(N_points); % batch size should be modified according to available GPU
end
batch_size = min(batch_size, prod(N_points));

% for storing sse
for i = 1:numel(grid_points)
    for j = 1:i
        if j==i
            SSE{i,j} = inf*ones(N_points(i),1);
        else
            SSE{i,j} = inf*ones(N_points(i),N_points(j));
        end
    end
end

% process batches sequentially
fig_handle = figure;
sse_best = inf;
for k = 1:ceil(prod(N_points)/batch_size)
    % get indices for each parameter for current batch
    inds = 1+(k-1)*batch_size:min(k*batch_size,prod(N_points));
    subs = cell(numel(grid_points),1);
    [subs{:}] = ind2sub(N_points,inds);

    % get associated parameter values
    pars = zeros(numel(grid_points),length(inds));
    for n = 1:numel(grid_points)
        pars(n,:) = grid_points{n}(subs{n});
    end

    % evaluate signals
    sse = zeros(1,size(pars,2));
    parfor i = 1:size(pars,2)
        sse(i) = sum((S_fun(pars(:,i)) - y).^2);
    end
  
 

    % store best encountered parameter combination
    [sse_new,best_index] = min(sse);
    if sse_new<sse_best
        sse_best = sse_new;
        pars_best = pars(:,best_index);
    end

    % store sse in bins
    for i = 1:numel(grid_points)
        for j = 1:i
            if j==i
                sorted_rows = sortrows(cat(1,subs{i},sse).');
                [min_subs,min_inds] = unique(sorted_rows(:,1),'rows');
                for n = 1:size(min_subs)
                    SSE{i,j}(min_subs(n,1)) = min(sorted_rows(min_inds(n),2),SSE{i,j}(min_subs(n,1)));
                end
            else
                sorted_rows = sortrows(cat(1,subs{i},subs{j},sse).');
                [min_subs,min_inds] = unique(sorted_rows(:,1:2),'rows');
                for n = 1:size(min_subs)
                    SSE{i,j}(min_subs(n,1),min_subs(n,2)) = min(sorted_rows(min_inds(n),3),SSE{i,j}(min_subs(n,1),min_subs(n,2)));
                end
            end
        end
    end

    % plot
    close(fig_handle)
    fig_handle = figure;
    colormap(flip(colormap))
    for i = 1:numel(grid_points)
        for j = 1:numel(grid_points)
            axes()
            if i~=j
                if i>j
                    map = SSE{i,j};
                else
                    map = SSE{j,i};
                    map = map';
                end
                map = length(y)*log(map/sse_best);
                s = pcolor(grid_points{j},grid_points{i},map);
                s.EdgeAlpha = 0;
                set(gca,'CLim',[0 10])
            else
                map = SSE{i,j};
                map = length(y)*log(map/sse_best);
                plot(grid_points{i},map(:),'o-')
                grid on
            end
            hold on, box on
            if i~=j
                plot(pars_best(j),pars_best(i),'r.','MarkerSize',10)
            else
                xline(pars_best(i),'k--');
            end
            if j~=1 && i~=j
                set(gca,'YTickLabels',[])
            end
            if i~=numel(grid_points)
                set(gca,'XTickLabels',[])
            end
            if j==size(pars,1) && i==1
                h = colorbar;
                ylabel(h,"$\Delta$IC")
            end
            if j==1
                ylabel(labels(i))
            end
            if i==numel(grid_points)
                xlabel(labels(j))
            end
        end
    end
    fig = myFigure(gcf,'size',[16 16],'pos',[1 1]);
    fig.subPlots([i i],'aspect_ratio',1,...
        'x1',1.3,'x2',0.6,'x3',2.0,...
        'y1',1.2,'y2',0.4,'y3',0.3)
    h.TickLabels{end} = "$\geq\,$" + h.TickLabels{end};
    pause(0.1)
end

