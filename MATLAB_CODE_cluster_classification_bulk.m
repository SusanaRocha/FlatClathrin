%% set file path: manually select folder. Load data.
fp = uigetdir;
cd(fp)
list=dir('*_DBscan.mat'); 
list={list.name}; 
 
th = 0.5035; % global threshold of the classification model. Set based on 1247 classical and 1042 alternative CCSs that were manually selected, ensuring that maximally 10% of the classical CCSs are misclassified    
flg = 1;

%% analyse files on by one
celldata = table('Size', [size(list, 2), 24], 'VariableTypes',{'string','single','single','single','single','single','single','single','single','single','single', 'single', 'single','single','single','single','single','single','single','single','single','single','single','string'}, ...
    'VariableNames', {'filename','nrFCL_nrclusters','nrCP_nrclusters','nrFCL_cell','nrCP_cell','nrFCL_nrCP','areaFCL_clathrinarea','areaCP_clathrinarea','areaFCL_cell','areaCP_cell','areaFCL_areaCP', 'nrptsFCL_totalnrpts','nrptsCP_totalnrpts','nrptsFCL_cell','nrptsCP_cell','nrptsFCL_nrptsCP','totalnrcluster','totalnrcluster_cellarea','cellarea','clathrinarea','clathrinarea_cellarea','totalnrpts','totalnrpts_cellarea','category'});
 
for f=1:size(list, 2)
    fn = list{1,f};
    file = [fp filesep fn];
    [prop_cl, cellarea] = get_prop_clusters(file,1, th); % runs separate MATLAB_CODE_get_prop_clusters, which calculates cluster properties, classifies each cluster as pit or lattice, and plots the result in an image
    if flg ==1
       save([file(1:strfind([fp fn], '.mat')) '_th_' num2str(th) '_individual_clusters_classified.mat'], 'prop_cl', 'cellarea');
    end
    ind_CP = strcmp(prop_cl.type(:), 'pit');
    ind_FCL = strcmp(prop_cl.type(:), 'lattice');
    celldata.filename(f) = fn;
   
    % per cell: number of clusters
    celldata.nrFCL_nrclusters(f) = sum(ind_FCL)./size(prop_cl, 1);
    celldata.nrCP_nrclusters(f) = sum(ind_CP)./size(prop_cl, 1);
    celldata.nrFCL_cell(f) = sum(ind_FCL)./cellarea;
    celldata.nrCP_cell(f) = sum(ind_CP)./cellarea;
    celldata.nrFCL_nrCP(f) = sum(ind_FCL)./sum(ind_CP);
  
    % per cell: area of clusters
    celldata.areaFCL_clathrinarea(f) = sum(prop_cl.area(ind_FCL))./sum(prop_cl.area(:));
    celldata.areaCP_clathrinarea(f) = sum(prop_cl.area(ind_CP))./sum(prop_cl.area(:));
    celldata.areaFCL_cell(f) = sum(prop_cl.area(ind_FCL))./cellarea;
    celldata.areaCP_cell(f) = sum(prop_cl.area(ind_CP))./cellarea;
    celldata.areaFCL_areaCP(f) = sum(prop_cl.area(ind_FCL))./sum(prop_cl.area(ind_CP));
   
    % per cell: nrpts of clusters
    celldata.nrptsFCL_totalnrpts(f) = sum(prop_cl.nrpts(ind_FCL))./sum(prop_cl.nrpts(:));
    celldata.nrptsCP_totalnrpts(f) = sum(prop_cl.nrpts(ind_CP))./sum(prop_cl.nrpts(:));
    celldata.nrptsFCL_cell(f) = sum(prop_cl.nrpts(ind_FCL))./cellarea;
    celldata.nrptsCP_cell(f) = sum(prop_cl.nrpts(ind_CP))./cellarea;
    celldata.nrptsFCL_nrptsCP(f) = sum(prop_cl.nrpts(ind_FCL))./sum(prop_cl.nrpts(ind_CP));
    
    % per cell: additional data
    celldata.totalnrcluster(f) = size(prop_cl, 1);
    celldata.totalnrcluster_cellarea(f) = size(prop_cl, 1)./cellarea;

    celldata.cellarea(f) = cellarea;
    celldata.clathrinarea(f) =  sum(prop_cl.area(:));
    celldata.clathrinarea_cellarea(f) = sum(prop_cl.area(:))./cellarea;

    celldata.totalnrpts(f) = sum(prop_cl.nrpts(:));
    celldata.totalnrpts_cellarea(f) = sum(prop_cl.nrpts(:))./cellarea;
end
    
save('all_files_cell_classification.mat','celldata') 
