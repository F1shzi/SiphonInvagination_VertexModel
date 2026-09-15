clc; clear; close all;

tic; disp(fix(clock)); 


hpf = [13 16 18];   const.hpf  = hpf;   hpf_Equil = 20;
Tamp = 100000;      const.Tamp = Tamp;
EquilStep = (hpf(end) - hpf(1)) * Tamp;
const.EquilStep = EquilStep;

DataStep  = 20000;  const.DataStep  = DataStep;
const.dt = 0.001;  


%%  ************** Model parameter *************** %
const.num = 50;   num = const.num;
const.area0 = 1;
const.W0 = 0.8;
const.H0 = const.area0 / const.W0;
const.R0 = const.num*const.W0/(2*pi);

const.Ka = 1;             % cell area modulus
const.Kc = 1e-3;          % passive contractility
const.Kb = 1e-4;          % tissue bending modulus
const.Km = 1.6e-3;        % active tension coefficient
const.Km_Lateral = 0.12;  % active tension coefficient
const.KaL = 1e-3;         % Tissue area modulus
const.areaL = pi*(const.R0-const.H0/2)^2;

const.GammaA = ones(num,1);   
const.GammaB = ones(num,1);   
const.GammaL = ones(num,1);

Elems = Mesh(const);
Elems.Myosin_jj = [1 1 1];
Force = ForceGroup(Elems, const);
VI = Elems.VI + const.dt * Force;
Elems.VI = VI; Elems = CellSize(Elems, const);
const.LengthScale = 11/(max(Elems.VI(Elems.Face_v0(Elems.CellIndex==0,:),2))-min(Elems.VI(Elems.Face_v0(Elems.CellIndex==0,:),2)));  
LengthScale = const.LengthScale;

Myosin0 = 1;   % Experimental Myosin parameters
MyosinB = Myosin0 .* [1; 1;     1;     1;     1;     1;     1;     1;     1;     1;     1    ];
MyosinA = MyosinB .* [1; 1.771; 3.148; 4.171; 4.909; 4.600; 3.039; 3.611; 4.087; 3.890; 4.472];
MyosinL = MyosinB .* [1; 1.555; 2.821; 2.578; 2.115; 3.785; 5.268; 6.560; 5.710; 5.013; 4.892];

Tlist_0 = Tamp * ((hpf(1) : 0.5 : hpf(end)) - hpf(1));
Tlist_1 = 1 : 1 : EquilStep;
MyosinA_T = interp1(Tlist_0, MyosinA, Tlist_1, 'linear')';
MyosinL_T = interp1(Tlist_0, MyosinL, Tlist_1, 'linear')';
MyosinB_T = interp1(Tlist_0, MyosinB, Tlist_1, 'linear')';
Elems.Myosin_t = [MyosinA_T ,MyosinL_T, MyosinB_T];

sigma_cell = 1.8;  % 1.5
w_gauss = exp(-(Elems.CellIndex.^2) / (2*sigma_cell^2));
Elems.Myosin_c = ones(num,3);

GammaL_index = [zeros(3*Tamp,1); linspace(0,1,1.5*Tamp)'; ones(floor(0.5*Tamp),1)]; % Lateral Tension increases


% Experimental measurements
Exp_Hpf = (13.5 : 0.5 : 18)';
Exp_InvagDepth = [0; 0; 0; 0; 0.05; 0.751619048; 1.807636364; 2.8173; 4.238761905; 6.5409];
Exp_CentHeight = [10.92961538; 12.58225; 13.4649; 14.8203; 14.94710526; 14.83519048; 14.17113636; 12.8881; 11.70142857; 11.58105];


%% Data Prepare
Data.Elemsstep  = cell(EquilStep / DataStep+1 ,1 );
Data.conststep  = cell(EquilStep / DataStep+1 ,1 );
Data.CenterInvagDepth = zeros(EquilStep,1);
Data.CenterCellHeight = zeros(EquilStep,1);
Data.CenterCellMyosin = zeros(EquilStep,3);

[~,VertexBetweenCell_34,~] = intersect(sort(Elems.NCI,2),sort([find(Elems.CellIndex==-3) find(Elems.CellIndex==-4) inf],2),'rows');
[~,VertexBetweenCell34,~] = intersect(sort(Elems.NCI,2),sort([find(Elems.CellIndex==3) find(Elems.CellIndex==4) inf],2),'rows');

%% Simulations
disp('#       Iteration     Time');
for jj = 1:EquilStep
    
    % Myosin update
    Myosin_jj = Elems.Myosin_t(jj,:);   % 1x3
    Myosin_Cell_jj = Myosin0 .* ones(num,3) + w_gauss * (Myosin_jj - Myosin0.*ones(1,3));
    Elems.Myosin_c = Myosin_Cell_jj;
    const.GammaA = const.Km .* (Myosin_Cell_jj(:,1)-Myosin0);  % const.GammaA = const.Km .* Myosin_Cell_jj(:,1);
    const.GammaB = const.Km .* (Myosin_Cell_jj(:,3)-Myosin0);  % const.GammaB = const.Km .* Myosin_Cell_jj(:,3);
    const.GammaL = const.Km_Lateral*GammaL_index(jj)*const.Km  .* (Myosin_Cell_jj(:,2)-Myosin0);  % 0.12 for ab
    
    % Vertex motion
    Force = ForceGroup(Elems, const);
    VI = Elems.VI + const.dt * Force;
    Elems.VI = VI;
    Elems = CellSize(Elems, const);
    
    % Save Data
  	if mod(jj,DataStep) == 0

        t0 = toc/60;  
        disp(['#        ',num2str(jj),'      ',num2str(t0,'%5.2f'),' min.']);

        Data = SaveData(Data,jj/DataStep,Elems,const);

    end

    Data.CenterCellHeight(jj)   = max(Elems.VI(Elems.Face_v0(Elems.CellIndex==0,:),2))-min(Elems.VI(Elems.Face_v0(Elems.CellIndex==0,:),2));
    Data.CenterCellMyosin(jj,:) = Elems.Myosin_c(Elems.CellIndex==0,:);
    Data.CenterInvagDepth(jj)   = max([(Elems.VI(VertexBetweenCell_34,2)+Elems.VI(VertexBetweenCell34,2))/2-max(Elems.VI(Elems.Face_v0(Elems.CellIndex==0,:),2)),0]);
       
end

%% Final shape after enough time
for jj = EquilStep+1 : EquilStep+hpf_Equil*Tamp

    Force = ForceGroup(Elems, const);
    VI = Elems.VI + const.dt * Force;
    Elems.VI = VI;
    Elems = CellSize(Elems, const);
    
  	if mod(jj,DataStep) == 0
        t0 = toc/60; 
        disp(['#        ',num2str(jj),'      ',num2str(t0,'%5.2f'),' min.']);
    end
    
end

Data = SaveData(Data,EquilStep / DataStep+1,Elems,const);

disp(fix(clock));  % time display



FileName = 'Results_control';
save([FileName,'.mat'],'Data','const')

