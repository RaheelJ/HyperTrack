function Done = Add_Paths()
    try
        addpath(genpath('Optimization/Optimal Control Tools/ICLOCS-master/src'))
		addpath(genpath('Optimization/Optimal Control Tools/OPTI-master'))
        addpath('Optimization')
        addpath('Configuration')
        Done = 1;
    catch
        Done = 0;
    end
end