function Elems = CellSize(Elems, const)

Nc = const.num;
VI = Elems.VI;
Face_v0 = Elems.Face_v0;
Nm = size(Face_v0,2);

% CALCULATE CELL AREA
verts = zeros(Nm,2*Nc);
verts(:,1:2:end) = reshape(VI(reshape(Face_v0',[],1),1),Nm,[]);
verts(:,2:2:end) = reshape(VI(reshape(Face_v0',[],1),2),Nm,[]);

vx = verts(:,1:2:end);   % Nm-by-Nc array
vy = verts(:,2:2:end);   % Nm-by-Nc array

area = polyarea(vx,vy)';

% CALCULATE CELL PERIMETER
vi = cat(2,reshape(vx,[],1),reshape(vy,[],1));
viR = cat(2,reshape(vx([2:end 1],:),[],1),reshape(vy([2:end 1],:),[],1));

vertsi_R = vi-viR;  
di_R = sqrt(sum(vertsi_R.^2,2));
perim = sum(reshape(di_R,Nm,[]),1)';   % Nc-by-1 array of perimeter of each cell

% Lumen area
vxL = VI(Elems.BasalVertices,1);
vyL = VI(Elems.BasalVertices,2);
areaL = polyarea(vxL,vyL)';

% Cell center
VI_all = VI(reshape(Face_v0',[],1),:);
Cgeo_x = mean(reshape(VI_all(:,1),4,[]),1);
Cgeo_y = mean(reshape(VI_all(:,2),4,[]),1);

%% 
Elems.Area = area;
Elems.Perim = perim;
Elems.AreaL = areaL;
Elems.CellCenter = [Cgeo_x',Cgeo_y'];

end