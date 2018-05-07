function OPG_Goats = Simu_OPG(OPG,nGoats)
%% Simulate OPG per animal
U = unique(OPG(:,2));
nt = randi(numel(U));
M = OPG(OPG(:,2) == U(nt),1);
M = M(randperm(length(M)));

OPG_Goats = M(1:nGoats);


