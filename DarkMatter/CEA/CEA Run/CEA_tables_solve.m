function [result,cstar,Tcc,Te,Pe,gamm_e,Isp,Ma_e,rho_e,R_e] = ...
    CEA_tables_solve(num_OF, OF, Pcc, SupAr, PcPe, fuel, ox)

    result = cell(40, num_OF);  % Cells to store data
    cstar = zeros(1,num_OF);
    Tcc = zeros(1,num_OF);
    Te = zeros(1,num_OF);
    Pe = zeros(1,num_OF);
    gamm_e = zeros(1,num_OF);
    Isp = zeros(1,num_OF);
    Ma_e = zeros(1,num_OF);
    rho_e = zeros(1,num_OF);    
    
    for i = 1:num_OF
        % Compute CEA solution
        ox{1, 3} = 278;
        [out] = CEA_Run(OF(i), Pcc, SupAr, PcPe, fuel, ox);

        % convert struct to cell for storage
        out_cell = struct2cell(out);

        % Add result to solution matrix
        result(1:length(out_cell),i) = out_cell;
        
        %name = fprintf('OF_%i',i);
        cstar(1,i)  = out.CSTAR.THROAT;
        Tcc(1,i)    = out.T.CHAMBER;
        Te(1,i)     = out.T.EXIT1;
        Pe(1,i)     = out.P.EXIT1 * 100000;
        gamm_e(1,i) = out.GAMMAs.EXIT1;
        Isp(1,i)    = out.Isp.EXIT1 / 9.81;
        Ma_e(1,i)   = out.MachNumber.EXIT1;
        rho_e(1,i)   = out.RHO.EXIT1;
    end
R_e = Pe ./ (rho_e .* Te);
end