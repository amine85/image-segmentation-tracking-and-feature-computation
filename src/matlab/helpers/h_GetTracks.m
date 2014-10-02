function [trackedImage] = h_GetTracks(labImage,nFrame)

% this modified tracking include inter-frames connections where frames are missing %
% global nFrame 

% count and initialize vertices loop %
% fprintf('entering count vertex loop\n');
vertex_list = cell(nFrame,1);
num_vertices = 0;

for tt = 1:nFrame
    currImg = uint16(labImage(:,:,tt));
    currIds = unique(currImg(:));
    currIds = currIds(currIds>0);
    currNumIds = length(currIds);
    
    curr_vert_list = cell(currNumIds,1);
    num_vertices = num_vertices + currNumIds;
    vertex_list{tt} = curr_vert_list;
end
% fprintf('toatl number of vertices is: %d \n',num_vertices);

%%%%%%%%%%%%%%%%%%% create vertex properties loop %%%%%%%%%%%%%%%%%%
% fprintf('creating vertex list\n');


tic;
vert_count = 1;
curr_frame = 1;
while(curr_frame<=nFrame)
    currLabImg = uint16(labImage(:,:,curr_frame));
    currIds = unique(currLabImg(:));
    currIds = currIds(currIds>0);
    if(isempty(currIds));
        curr_frame = curr_frame+1;
        continue;
    end    
%     curr_frame
    for ii = 1:length(currIds)
            curr_id = currIds(ii);
            binImg = (currLabImg == curr_id);            
            props = regionprops(binImg, 'Centroid','Area','PixelList','PixelIdxList');  

            % vertex properties %
            vertex.hull_id =  curr_id;
            vertex.hull_centroid = props.Centroid;
            vertex.hull_volume = props.Area;
            vertex.hull_pixList = props.PixelList;
            vertex.hull_linPixList = props.PixelIdxList;
            vertex.time = curr_frame;            
            vertex.vert_id = vert_count;
            vertex.new_id = -1;
            vertex.visited = false;
            vertex.in_edge_indx_list = [];
            vertex.out_edge_indx_list = [];            
            vertex_list{curr_frame}{ii} = vertex;
            vert_count = vert_count + 1;
    end   
    curr_frame = curr_frame+1;
    
end
op_time = toc;
% fprintf('vertex list time: %d \n',op_time);
            
%%%%%%%%%%%%%%%%%%% estimate cost paramters %%%%%%%  

% % % curr_frame = 1;
% % % finished = false;
% % % while(curr_frame<nFrame)
% % %     M = length(vertex_list{curr_frame});
% % %     if(M==0)
% % %         curr_frame = curr_frame+1;
% % %         continue;
% % %     end
% % %     
% % %     next_frame = curr_frame+1;
% % %     
% % %     while(isempty(vertex_list{next_frame}))
% % %         if(next_frame == nFrame)
% % %             finished = true;
% % %             break;
% % %         end
% % %         next_frame = next_frame+1;
% % %     end
% % %     
% % %     if(finished)
% % %         break;
% % %     end
% % %     
% % %     N = length(vertex_list{next_frame});
% % % 
% % %     currCost = zeros(M*N,1);
% % %     counter = 1;
% % %     for ii = 1:M
% % %         for jj = 1:N
% % %            cen_dist = GetCentDistance(vertex_list{curr_frame}{ii}.hull_centroid,vertex_list{next_frame}{jj}.hull_centroid);
% % %            vol_dist = GetVolDistance(vertex_list{curr_frame}{ii}.hull_volume,vertex_list{next_frame}{jj}.hull_volume);
% % %            cc_dist = GetCCDistance(vertex_list{curr_frame}{ii}.hull_linPixList, vertex_list{next_frame}{jj}.hull_linPixList,...
% % %                                     vertex_list{curr_frame}{ii}.hull_pixList, vertex_list{next_frame}{jj}.hull_pixList);
% % %            cost = (10* cen_dist + 100* vol_dist + 1000* cc_dist)^2;
% % %            currCost(counter) = cost; 
% % %            counter = counter + 1;      
% % %         end
% % %     end
% % %     costList{curr_frame} = currCost;
% % %     curr_frame = next_frame;
% % % end



%%%%%%%%%%%%%%%%%%%%%%% add inter-frame edges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_edge_id = 1;                      % global edge id

curr_frame = 1;
finished = false;
edge_list = cell(nFrame-1,1);
while(curr_frame<nFrame)
    M = length(vertex_list{curr_frame});
    if(M==0)
%         edge_list{curr_frame} = [];
        curr_frame = curr_frame+1;
        continue;
    end
    
    next_frame = curr_frame+1;
    
    while(isempty(vertex_list{next_frame}))
        if(next_frame == nFrame)
            finished = true;
            break;
        end
        next_frame = next_frame+1;
    end
    
    if(finished)
        break;
    end
    
    N = length(vertex_list{next_frame});
    
    % get max cost %
%     maxCost = max(costList{curr_frame});
    edge_count = 1;
    curr_edge_list = cell(M*N,1);
    for ii = 1:M
        for jj = 1:N
        
           edge.in_vert_time = curr_frame;
           edge.in_vert_indx = ii;
           edge.out_vert_time = next_frame;
           edge.out_vert_indx = jj;
           edge.selected = false;
           edge.id = g_edge_id;
           
           % compute edge cost %
           cen_dist = GetCentDistance(vertex_list{curr_frame}{ii}.hull_centroid,vertex_list{next_frame}{jj}.hull_centroid);
           vol_dist = GetVolDistance(vertex_list{curr_frame}{ii}.hull_volume,vertex_list{next_frame}{jj}.hull_volume);
           cc_dist = GetCCDistance(vertex_list{curr_frame}{ii}.hull_linPixList, vertex_list{next_frame}{jj}.hull_linPixList,...
                                    vertex_list{curr_frame}{ii}.hull_pixList, vertex_list{next_frame}{jj}.hull_pixList);

%            cost = (10* cen_dist + 100* vol_dist + 1000* cc_dist)^2/maxCost;
           cost = 10* cen_dist + 100* vol_dist + 1000* cc_dist;
%            edge.cost = exp(-6*cost);
           edge.cost = cost;
           curr_edge_list{edge_count} = edge;
           edge_count = edge_count + 1;
           
            vertex_list{curr_frame}{ii}.out_edge_indx_list = [vertex_list{curr_frame}{ii}.out_edge_indx_list;g_edge_id];
            vertex_list{next_frame}{jj}.in_edge_indx_list = [vertex_list{next_frame}{jj}.in_edge_indx_list;g_edge_id];
            
           g_edge_id = g_edge_id + 1;

        end       
    end
    
    % normalize the cost from 0-1
    costs = zeros(1,length(curr_edge_list));
    for kk = 1:length(curr_edge_list)
        costs(kk) = curr_edge_list{kk}.cost;
    end
    maxCost = max(costs);

    
    for kk = 1:length(curr_edge_list)
        dum = (curr_edge_list{kk}.cost/(maxCost+eps))^2;
        curr_edge_list{kk}.cost =  exp(-6*dum);
    end

    edge_list{curr_frame} = curr_edge_list;  
    
    curr_frame = next_frame;

    
end% end of frame loop 
oper_time = toc;
% fprintf('constraint list time : %d \n',oper_time);
% fprintf('total number of edges : %d \n',g_edge_id);

if(g_edge_id==1)
    trackedImage = labImage;
    return;
end


% count the number of edges %
total_num_edges = 0;
% length(edge_list)
for tt = 1:nFrame - 1    
    total_num_edges = total_num_edges + length(edge_list{tt});
end
% fprintf('total number of edges is : %d \n',total_num_edges);


%%%%%%%%%%%%%%%%%%%%% generate cost function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
costFunction = zeros(total_num_edges,1);
for tt = 1:nFrame - 1    
    for ee = 1:length(edge_list{tt})
       edge_id =  edge_list{tt}{ee}.id;
       costFunction(edge_id) = edge_list{tt}{ee}.cost;        
    end
end

%%%%%%%%%%%%%%%%%%%%%% generate constraint matrix %%%%%%%%%%%%%%%%%%%%%%%%%%

% count the number of out_edges and in_edges %
num_out_edges_constraints = 0;
num_in_edges_constraints = 0;

for tt = 1:nFrame
    for vv = 1:length(vertex_list{tt})
        if(~isempty(vertex_list{tt}{vv}.out_edge_indx_list))
            num_out_edges_constraints = num_out_edges_constraints +1;
        end
        if(~isempty(vertex_list{tt}{vv}.in_edge_indx_list))
            num_in_edges_constraints = num_in_edges_constraints +1;
        end
    end
end
% fprintf('total number of out_edges constraints is : %d \n',num_out_edges_constraints);
% fprintf('total number of in_edges constraints is : %d \n',num_in_edges_constraints);

constraintMatrix = zeros(num_out_edges_constraints+num_in_edges_constraints,total_num_edges);

% add out_edge constraints %
constr_count = 1;
for tt = 1:nFrame
    for vv = 1:length(vertex_list{tt})
        if(~isempty(vertex_list{tt}{vv}.out_edge_indx_list))
            out_edge_indx = vertex_list{tt}{vv}.out_edge_indx_list;
            constraintMatrix(constr_count,out_edge_indx)=1;
            constr_count = constr_count + 1;
        end
    end
end

% add in_edge constraints %
for tt = 1:nFrame
    for vv = 1:length(vertex_list{tt})
        if(~isempty(vertex_list{tt}{vv}.in_edge_indx_list))
            in_edge_indx = vertex_list{tt}{vv}.in_edge_indx_list;
            constraintMatrix(constr_count,in_edge_indx)=1;
            constr_count = constr_count + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%% solve the integer program %%%%%%%%%%%%%%%%%%%%%%%%%%
 b = ones(num_out_edges_constraints+num_in_edges_constraints,1);
 options = optimset('Display','off');
 warning('off','optim:bintprog:NotSupported')
 [x,~,exitflag] = bintprog(-costFunction,constraintMatrix,b,[],[],[],options);
%  fprintf('integer program exited with: %d\n',exitflag);
 %  x = bintprog(costFunction,constraintMatrix,b);

 % label selected edges %
 for tt = 1:nFrame - 1    
     if(isempty( edge_list{tt}))
        continue;
     end

     for ee = 1:length(edge_list{tt})
        if(x(edge_list{tt}{ee}.id))
            edge_list{tt}{ee}.selected = true;            
        end
     end
 end
 
 %%%%%%%%%%%%%%%%%%%% relabel vertex ids %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('entering relabeling loop\n');
new_id = 1;
tic;

for tt = 1:nFrame - 1
    if(isempty( edge_list{tt}))
        continue;
    end
    
    for ee = 1:length(edge_list{tt})
        if (edge_list{tt}{ee}.selected)
            
            in_vert_time = edge_list{tt}{ee}.in_vert_time;
            in_vert_indx = edge_list{tt}{ee}.in_vert_indx;
            
            out_vert_time = edge_list{tt}{ee}.out_vert_time;
            out_vert_indx = edge_list{tt}{ee}.out_vert_indx;
            
            if(~vertex_list{in_vert_time}{in_vert_indx}.visited)
                vertex_list{in_vert_time}{in_vert_indx}.new_id = new_id;
                vertex_list{out_vert_time}{out_vert_indx}.new_id = new_id;
                vertex_list{in_vert_time}{in_vert_indx}.visited = true;    
                vertex_list{out_vert_time}{out_vert_indx}.visited = true;   
                new_id = new_id +1;
            else
                vertex_list{out_vert_time}{out_vert_indx}.new_id = vertex_list{in_vert_time}{in_vert_indx}.new_id;
                vertex_list{out_vert_time}{out_vert_indx}.visited = true;   
            end
        end
    end
end
process_time = toc;
% fprintf('relabeling time : %d \n',process_time);


 
% create tracked  images %
 new_id = new_id + 1;
trackedImage = zeros(size(labImage));
for tt = 1:nFrame
    currLabImg = uint16(labImage(:,:,tt));
    currTrackImg = zeros(size(currLabImg));
    for vv = 1:length(vertex_list{tt})
        if(vertex_list{tt}{vv}.visited)
             currTrackImg(currLabImg == vertex_list{tt}{vv}.hull_id) = vertex_list{tt}{vv}.new_id;
        else
             currTrackImg(currLabImg == vertex_list{tt}{vv}.hull_id) = new_id;
             new_id = new_id + 1;
        end
    end
    trackedImage(:,:,tt) = currTrackImg;
end



%%% utility functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % centroid distance %
    function dist = GetCentDistance(x,y)
        dist = sqrt(sum((x-y).^2));        
    end
    
    % volume distance %
    function voldist = GetVolDistance(v1,v2)
       vmax = max([v1,v2]);
       vmin = min([v1,v2]);
       voldist = (vmax-vmin)/vmax;        
    end

    % connected component distance %
    function cc_dist = GetCCDistance(linIdx1,linIdx2,coord1,coord2)
%         size(linIdx1)
        % compute overlap %
%         overlap = 0;
%         for i = 1:length(linIdx1)
%             for j = 1:length(linIdx2)
%                 if(linIdx1(i)==linIdx2(j))
%                    overlap = overlap +1;
%                 end               
%             end
%         end
        overlap = length(intersect( linIdx1, linIdx2 ));
%         if( overlap == overlap2 )
%             disp('yes')
%         end
        
        if(overlap>0)
           vmin = min([length(linIdx1),length(linIdx2)]);
           cc_dist = 1-(overlap/vmin);          
        else
            D = pdist2(coord1,coord2);
            cc_dist = min(D(:));            
        end
    end



end % end of function
































