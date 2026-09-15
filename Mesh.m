function Elems = Mesh(const)

num   = const.num;  

R = const.R0;
Ri = R - const.H0 - 0.2;
Ro = R + const.H0 + 0.2;

list = zeros(num,2);

theta = pi/2;
for i = 1 : num
    list(i,1) = R * cos(theta);
    list(i,2) = R * sin(theta);
    theta = theta + 2*pi / num;
end

Nm = round((5 * num)^ 0.5) ^ 2;
list_mi = zeros(Nm, 2);
list_mo = zeros(Nm, 2);
theta = 0;
for i = 1 : Nm
    list_mi(i,1) = Ri * cos(theta);
    list_mi(i,2) = Ri * sin(theta);
    list_mo(i,1) = Ro * cos(theta);
    list_mo(i,2) = Ro * sin(theta);
    theta = theta + 2*pi / Nm;
end

list = [list; list_mi; list_mo];

%% ******************** Set up  mesh ******************** %%
[vertex,index] = voronoin(list);

clist_v = index(1:num);

[sortvertex, ~] = sort(cat(2,clist_v{:}));
Frequencytable = tabulate(sortvertex);

vertex_nci3 = find(Frequencytable(:,2) == 3);
vertex_nci2 = find(Frequencytable(:,2) == 2);

vertex_new = [vertex_nci3 ; vertex_nci2];
VI = vertex(vertex_new,:);
Nv = size(VI,1);

vertex_nci3 = (1:1:size(vertex_nci3,1))';
vertex_nci2 = (size(vertex_nci3,1)+1:1:size(vertex_nci3,1)+size(vertex_nci2,1))';

%% ********************* Set up  NCI (Neighboring cell indices of each vertex)  ********************* %%
NCI = zeros(Nv,3);

vlist = cat(2,clist_v{:})';

cellidxs = cellfun(@(x,y) y.*ones(size(x)),clist_v,num2cell((1:length(clist_v))'),'uniformoutput',false);
cellidxs = cat(2,cellidxs{:})';   % get the cell index each vertex belongs to

for i = 1 : Nv
    vi = vertex_new(i);
    ncii = cellidxs(vlist==vi)';   % neighboring cell of each new vertex
                if vecnorm(VI(i,:),2,2) > R
                    NCI(i,:) = [ncii,10*num];
                else
                    NCI(i,:) = [ncii,-10*num];
                end
end    

cgeo = zeros(num,2);
Face_v = cell(num,1);

for c = 1 : num
    [~,~,vcnew] = intersect(clist_v{c}, vertex_new);
    Face_v{c} = vcnew;   % the index of vertices of ith cell
    cgeo(c,:) = sum(VI(vcnew,:),1) ./ size(vcnew,1);
end

%% ********************* Set up  NVI (Neighboring vertex indices of each vertex) ********************* %%
NCI0 = NCI;

nci_v = zeros(3*Nv,2);    % 2 neighbour cells define 1 neighbour vertex
nci_v(1:3:end,:) = NCI0(:,[1 2]);
nci_v(2:3:end,:) = NCI0(:,[2 3]);
nci_v(3:3:end,:) = NCI0(:,[3 1]);
nci_v = sort(nci_v,2);

dmat2 = pdist2(nci_v,nci_v);
dmat2(1:3*Nv+1:end) = Inf;
[~,list_nbor] = min(dmat2,[],2);
NVI = reshape(ceil(list_nbor/3),3,[])';

%% ****************** Reorder NCI, NVI in orders ****************** %
vigroup = zeros(4,2*Nv);   
vigroup(1,:) = reshape(VI',1,[]);
vigroup(2,:) = reshape(cgeo(NCI(:,1),:)',1,[]);
vigroup(3,:) = reshape(cgeo(NCI(:,2),:)',1,[]);

for v = 1:Nv
        vigroup(4,2*v-1:2*v) = 3 .* vigroup(1,2*v-1:2*v)...
            -  vigroup(2,2*v-1:2*v) -  vigroup(3,2*v-1:2*v);
end

rJ1_i = reshape(vigroup(2,:)-vigroup(1,:),2,[])';  
rJ2_i = reshape(vigroup(3,:)-vigroup(1,:),2,[])';  
rJ3_i = reshape(vigroup(4,:)-vigroup(1,:),2,[])';  

judge1 = Cross2D(rJ1_i,rJ2_i);
judge2 = Cross2D(rJ2_i,rJ3_i);
judge3 = Cross2D(rJ3_i,rJ1_i);

judge = cat(2,judge1<0,judge2<0,judge3<0);   % find the index of vertices whose neighbor cells are clockwise ordered
judge = sum(judge,2);
flag = judge > 1;

NCI(flag,1:3) = NCI(flag,[1 3 2]);
NVI(flag,1:3) = NVI(flag,[3 2 1]);


%% ************** Reorder the vertex order of cells ************** %
for c = 1: num
    VertsNotInOrder =  Face_v{c};
    Nc = size(VertsNotInOrder,1);
    VertsInOrder = VertsNotInOrder;
    VertsInOrder(2:Nc) = zeros(Nc-1,1);
    if Nc > 3
        nextv = intersect(NVI(VertsNotInOrder(1),:),VertsNotInOrder);
        k = 2;
        while k ~= Nc
            VertsInOrder(k) = nextv(1);
            nextv = intersect(NVI(VertsInOrder(k),:),VertsNotInOrder);
            nextv = setdiff(nextv,VertsInOrder);
            k = k + 1;
        end
        VertsInOrder(k) = nextv(1);
    end
    VI_VertsInOrder = VI(VertsInOrder,:);
    rI_j1 = VI_VertsInOrder(1,:) - cgeo(c,:);
    rI_j1 = rI_j1 ./ vecnorm(rI_j1,2,2);
    rI_j2 = VI_VertsInOrder(2,:) - cgeo(c,:);
    rI_j2 = rI_j2 ./ vecnorm(rI_j2,2,2);
    judge = Cross2D(rI_j1,rI_j2);
    if judge < 0
        VertsInOrder = flipud(VertsInOrder);
    end

    Face_v{c} = VertsInOrder';
    
end

%% *************** Identify lateral edges *************** %
NCI(NCI==10*num) = inf;  % Apical boundary
NCI(NCI==-10*num) = -1;  % Basal  boundary

nvi = reshape(NVI',[],1);
nvi0 = reshape(((1:1:Nv)'.*ones(Nv,3))',[],1);
neighboringcell = sort([NCI(nvi0,:) , NCI(nvi,:)], 2);
BasalEdge   = neighboringcell(:,2) < 0;
ApicalEdge  = neighboringcell(:,5)>=Inf;
LateralEdge0 = double(BasalEdge)+double(ApicalEdge)==0;
LateralEdge0 = reshape(LateralEdge0,3,[])';     % All lateral edge

LateralEdge = unique(sort([(1:1:Nv)',nvi(reshape(LateralEdge0',[],1))],2),'rows');

LateralEdge_NeighborCell = zeros(size(LateralEdge,1),2);
for e=1:size(LateralEdge,1)
    LateralEdge_NeighborCell(e,:) =  setdiff(intersect(NCI(LateralEdge(e,1),:),NCI(LateralEdge(e,2),:)),[-1,inf]);
end
LateralEdge = [LateralEdge,LateralEdge_NeighborCell];

%% *************** Identify apical/basal vertices *************** %
nci = sort(NCI,2);
Va = find(nci(:,3)==inf);   % Apical vertices
Vb = find(nci(:,1)==-1);    % Basal vertices
Nva = length(Va);           
Nvb = length(Vb);

List_Va = zeros(Nva,1); % connect apical vertices in order
List_Va(1) = Va(1);
q0 = intersect(NVI(Va(1),:),Va);
List_Va(2) = q0(1);
for i = 3 : length(Va)
    List_Va(i)  = setdiff(intersect(NVI(List_Va(i-1),:),Va),List_Va(i-2));
end

List_Vb = zeros(Nvb,1); % connect basal vertices in order
List_Vb(1) = Vb(1);
q0 = intersect(NVI(Vb(1),:),Vb);
List_Vb(2) = q0(1);
for i = 3 : length(Vb)
    List_Vb(i)  = setdiff(intersect(NVI(List_Vb(i-1),:),Vb),List_Vb(i-2));
end

ApicalEdge = [List_Va, [List_Va(2:end);List_Va(1)]];
ApicalEdge_NeighborCell = zeros(size(List_Va,1),1);
for e=1:size(List_Va,1)
    ApicalEdge_NeighborCell(e) =  setdiff(intersect(NCI(ApicalEdge(e,1),:),NCI(ApicalEdge(e,2),:)),[-1,inf]);
end
ApicalEdge = [ApicalEdge,ApicalEdge_NeighborCell];

BasalEdge  = [List_Vb, [List_Vb(2:end);List_Vb(1)]];
BasalEdge_NeighborCell = zeros(size(List_Vb,1),1);
for e=1:size(List_Vb,1)
    BasalEdge_NeighborCell(e) =  setdiff(intersect(NCI(BasalEdge(e,1),:),NCI(BasalEdge(e,2),:)),[-1,inf]);
end
BasalEdge = [BasalEdge,BasalEdge_NeighborCell];

%% Define Cell index
%   top-center cell: 0
%   left side:  -1, -2, -3, ...
%   right side: +1, +2, +3, ...

theta = atan2(cgeo(:,2), cgeo(:,1));     % angle in [-pi, pi]
phi = pi/2 - theta;

% wrap phi into (-pi, pi]
phi(phi <= -pi) = phi(phi <= -pi) + 2*pi;
phi(phi >   pi) = phi(phi >   pi) - 2*pi;

[~, idx0] = min(abs(phi)); % find the top-center cell: the one with minimal |phi|
c0 = idx0; 

left_cells  = setdiff(find(phi < 0),c0); % split cells into left and right according to sign of phi
right_cells = setdiff(find(phi > 0),c0);

[~, il] = sort(abs(phi(left_cells)), 'ascend'); % sort left side: from top to bottom along left arc
left_sorted = left_cells(il);

[~, ir] = sort(abs(phi(right_cells)), 'ascend'); % sort right side: from top to bottom along right arc
right_sorted = right_cells(ir);

% initialize CellIndex
CellIndex = nan(num,1);
CellIndex(c0) = 0;

for k = 1:length(left_sorted)          % assign negative indices on the left
    CellIndex(left_sorted(k)) = -k;
end

for k = 1:length(right_sorted)         % assign positive indices on the right
    CellIndex(right_sorted(k)) = k;
end

%% ********************** Save in Elements ********************** %
Facesize = cellfun(@(x) length(x),Face_v);
Nm = max(Facesize);
func = @(x) x([1:end ones(1,Nm-end)]);
Face_v0 = cellfun(func,Face_v,'UniformOutput',false);
Face_v0 = cell2mat(Face_v0);

Elems.VI = VI;
Elems.NCI = NCI;
Elems.NVI = NVI;
Elems.Face_v = Face_v;
Elems.Face_v0 = Face_v0;
Elems.CellCenter = cgeo;

Elems.ApicalVertices = List_Va;
Elems.BasalVertices = List_Vb;
Elems.LateralEdge = LateralEdge;
Elems.ApicalEdge = ApicalEdge;
Elems.BasalEdge = BasalEdge;

Elems.CellIndex   = CellIndex;

Elems = CellSize(Elems, const);



%% ************************** End ******************************* %

% SUBFUNCTIONS
    function c = Cross2D(a,b)
        
        c = a(:,1).*b(:,2)-a(:,2).*b(:,1);
        
    end



end