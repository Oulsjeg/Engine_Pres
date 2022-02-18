function addfxnpath()

    if ismac || isunix

            addpath(strcat([pwd, '/CEA']));
            addpath(strcat([pwd, '/CEA/CEA Run']));
            addpath(strcat([pwd, '/src']));

    elseif ispc

            addpath(strcat([pwd, '\CEA']));
            addpath(strcat([pwd, '\CEA\CEA Run']));
            addpath(strcat([pwd, '\src']));

    else
        disp('Platform not supported')
    end
end
