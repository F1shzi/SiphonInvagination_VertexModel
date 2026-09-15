function  Force = ForceGroup(Elems, const)


Ka  = const.Ka;   
Kc  = const.Kc;   
Kb  = const.Kb;   
KaL = const.Ka;   

GammaA = [const.GammaA ; 0];   
GammaB = [const.GammaB ; 0];   
GammaL = [const.GammaL ; 0];   

area0 = const.area0;   
areaL = const.areaL;

Nc = const.num;
VI = Elems.VI;
Nv = size(VI,1);
NCI = Elems.NCI;
NVI = Elems.NVI;

%% =================== Area elasticity & Lumen area ================== %
nci = reshape(NCI',[],1);
nvi0 = reshape(((1:1:Nv)'.*ones(Nv,3))',[],1);
NVI1 = NVI(:,[3 1 2]);       nvi1 = reshape(NVI1',[],1);
NVI2 = NVI(:,[1 2 3]);       nvi2 = reshape(NVI2',[],1);
VI3 = [VI , zeros(Nv,1)];

rj1_j2 = VI3(nvi2,:) - VI3(nvi1,:);

rj2_i = VI3(nvi0,:) - VI3(nvi2,:);
ri_j1 = VI3(nvi1,:) - VI3(nvi0,:);
k = cross(rj2_i,ri_j1);                              
k = k ./ vecnorm(k,2,2);

Area = [Elems.Area   ; Elems.AreaL; 0 ];
area = [repmat(area0,Nc,1);  areaL; 0 ];
Ka_all = [repmat(Ka,Nc,1);     KaL; 0 ];

nci(nci == -1 ) = Nc + 1;
nci(nci == Inf) = Nc + 2;

Coefficient0_Fa = 0.5 .* ( - Ka_all .*( Area - area) );
Coefficient_Fa  = Coefficient0_Fa(nci);

Fa = 2 .* Coefficient_Fa .* cross(k , rj1_j2);
Fa = Fa(1:3:end,:) +  Fa(2:3:end,:) +  Fa(3:3:end,:);
Fa = Fa(:,1:2);
 
%% =================== Cortical contraction ================== %
rj1_i = VI(nvi0,:) - VI(nvi1,:);  rj1_i = rj1_i ./ vecnorm(rj1_i,2,2);
rj2_i = VI(nvi0,:) - VI(nvi2,:);  rj2_i = rj2_i ./ vecnorm(rj2_i,2,2);

Coefficient0_Fc = - 0.5 .* Kc .* [ Elems.Perim ; 0 ; 0];
Coefficient_Fc = Coefficient0_Fc(nci);

Fc = 2 .* Coefficient_Fc .* (rj1_i + rj2_i);
Fc = Fc(1:3:end,:) +  Fc(2:3:end,:) +  Fc(3:3:end,:);

%% =================== Active Tension energy ================== %
nvi = reshape(NVI',[],1);
nvi0 = reshape(((1:1:Nv)'.*ones(Nv,3))',[],1);
neighboringcell = sort([NCI(nvi0,:) , NCI(nvi,:)], 2);
NCI1 = NCI(:,[1 2 3]);       nci1 = reshape(NCI1',[],1);
NCI2 = NCI(:,[2 3 1]);       nci2 = reshape(NCI2',[],1);

BasalEdge   = neighboringcell(:,2) < 0;
ApicalEdge  = neighboringcell(:,5)>=Inf;
LateralEdge = double(BasalEdge)+double(ApicalEdge)==0;  

ri_j = VI(nvi0,:)-VI(nvi,:);
ri_j_L = vecnorm(ri_j,2,2);
ri_j = ri_j ./ ri_j_L;

nci1(nci1 == -1 ) = Nc + 1; nci1(nci1 == Inf) = Nc + 1;
nci2(nci2 == -1 ) = Nc + 1; nci2(nci2 == Inf) = Nc + 1;

Coefficient_Ft = zeros(size(nvi,1),1);
Coefficient_Ft(LateralEdge) = - 0.5.*(GammaL(nci1(LateralEdge))+GammaL(nci2(LateralEdge)));
Coefficient_Ft(ApicalEdge)  = - (GammaA(nci1(ApicalEdge))+GammaA(nci2(ApicalEdge)));
Coefficient_Ft(BasalEdge)   = - (GammaB(nci1(BasalEdge))+GammaB(nci2(BasalEdge)));


Ft = 2 .* ri_j_L .* Coefficient_Ft .* ri_j;
Ft = Ft(1:3:end,:) +  Ft(2:3:end,:) +  Ft(3:3:end,:);


%% =================== Bending energy ================== %
% apical bending force
List_Va = Elems.ApicalVertices;
List_Va1 = [List_Va(2:end);List_Va(1)] ;            % index i+1
List_Va2 = [List_Va(3:end);List_Va(1:2)] ;          % index i+2
List_Va_1 = [List_Va(end);List_Va(1:end-1)] ;       % index i-1
List_Va_2 = [List_Va(end-1:end);List_Va(1:end-2)] ; % index i-2

ri   = VI(List_Va,:);       % ri
ri1  = VI(List_Va1,:);      % ri+1
ri2  = VI(List_Va2,:);      % ri+2
ri_1 = VI(List_Va_1,:);     % ri-1
ri_2 = VI(List_Va_2,:);     % ri-2

CosThetai = dot((ri1-ri),(ri-ri_1),2) ./ (vecnorm(ri1-ri,2,2) .* vecnorm(ri-ri_1,2,2));
CosThetai1 =  [CosThetai(2:end);CosThetai(1)] ;           % cos(theta)i+1
CosThetai_1 = [CosThetai(end);CosThetai(1:end-1)] ;       % cos(theta)i-1

Fb_Va = (CosThetai  -1).*((ri1+ri_1-2.*ri)./(vecnorm(ri1-ri,2,2).*vecnorm(ri_1-ri,2,2))+...
                           CosThetai.*((ri1-ri)./(vecnorm(ri1-ri,2,2).^2)+(ri_1-ri)./(vecnorm(ri_1-ri,2,2).^2)))+...
        (CosThetai_1-1).*((ri_1-ri_2)./(vecnorm(ri-ri_1,2,2).*vecnorm(ri_2-ri_1,2,2))+...
                           CosThetai_1.*((ri_1-ri)./(vecnorm(ri_1-ri,2,2).^2)))+...
        (CosThetai1 -1).*((ri1-ri2)./(vecnorm(ri-ri1,2,2).*vecnorm(ri2-ri1,2,2))+...
                           CosThetai1.*((ri1-ri)./(vecnorm(ri1-ri,2,2).^2)));
                       
% basal bending force
List_Vb = Elems.BasalVertices;
List_Vb1 = [List_Vb(2:end);List_Vb(1)] ;            % index i+1
List_Vb2 = [List_Vb(3:end);List_Vb(1:2)] ;          % index i+2
List_Vb_1 = [List_Vb(end);List_Vb(1:end-1)] ;       % index i-1
List_Vb_2 = [List_Vb(end-1:end);List_Vb(1:end-2)] ; % index i-2

ri   = VI(List_Vb,:);       % ri
ri1  = VI(List_Vb1,:);      % ri+1
ri2  = VI(List_Vb2,:);      % ri+2
ri_1 = VI(List_Vb_1,:);     % ri-1
ri_2 = VI(List_Vb_2,:);     % ri-2

CosThetai = dot((ri1-ri),(ri-ri_1),2) ./ (vecnorm(ri1-ri,2,2) .* vecnorm(ri_1-ri,2,2));
CosThetai1 =  [CosThetai(2:end);CosThetai(1)] ;          
CosThetai_1 = [CosThetai(end);CosThetai(1:end-1)] ;      

Fb_Vb = (CosThetai  -1).*((ri1+ri_1-2.*ri)./(vecnorm(ri1-ri,2,2).*vecnorm(ri_1-ri,2,2))+...
                           CosThetai.*((ri1-ri)./(vecnorm(ri1-ri,2,2).^2)+(ri_1-ri)./(vecnorm(ri_1-ri,2,2).^2)))+...
        (CosThetai_1-1).*((ri_1-ri_2)./(vecnorm(ri-ri_1,2,2).*vecnorm(ri_2-ri_1,2,2))+...
                           CosThetai_1.*((ri_1-ri)./(vecnorm(ri_1-ri,2,2).^2)))+...
        (CosThetai1 -1).*((ri1-ri2)./(vecnorm(ri-ri1,2,2).*vecnorm(ri2-ri1,2,2))+...
                           CosThetai1.*((ri1-ri)./(vecnorm(ri1-ri,2,2).^2)));   
                             
Fb = zeros(Nv,2);
Fb(List_Va,:) = - 2 * Kb .* Fb_Va;
Fb(List_Vb,:) = - 2 * Kb .* Fb_Vb;



%% *************************** End **************************** %

Force =  Fa + Fc + Fb + Ft ;


end