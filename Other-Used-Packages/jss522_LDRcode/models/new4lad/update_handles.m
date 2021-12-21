function [Fhandle,dFhandle] = update_handles(data_parameters)
Fhandle = F(@F4newlad,data_parameters);
dFhandle = dF(@dF4newlad_analitica,data_parameters);
