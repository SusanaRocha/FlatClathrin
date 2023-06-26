# FlatClathrin
<h2> READ ME </h2>
This document contains MATLAB code for <br>
•	Localizing single molecules in PALM-TIRF experiment in bulk <br>
•	Identifying  clusters and calculating certain cluster properties <br>
•	Classifying clusters as a pit or FCL <br>
<br>
<h2> MATLAB_CODE_Localizations_bulk </h2>
<h3>Input </h3><br>
•	.HIS files from single-molecule PALM-TIRF experiment. For one cell, two movies of each 5000 frames are imaged. The two movies have to have the same core name, followed by _mov1 or _mov2 for the analysis to work. Example: 20211126_KM12L4a_mEos3.2_03_mov1.HIS and 20211126_KM12L4a_mEos3.2_03_mov2.HIS (core name = 20211126_KM12L4a_mEos3.2_03) - <i>data can be supplied upon request </i> <br>
•	Localizer function (Dedecker et al. doi:10.1117/1.JBO.17.12. For download and installation, see: https://bitbucket.org/pdedecker/localizer/src/master/ ) <br>
<h3>Important code aspects </h3><br>
•	Bulk analysis for .HIS files <br>
•	First individual .HIS files are analyzed: <br>
o	Localizations are retrieved by using the Localizer function that fits a 2D Gaussian with PSF standard deviation factor 1.8 and intensity selection sigma factor 25. Localizations are saved in the variable pts: a table with 12 columns (see explanation Localizer function) where each row represents 1 localization. <br>
o	Localizations are plotted with a scale bar of 5 µm. LUT scale of these images can be adjusted in the script.<br>
•	Then, two movies with the same core name (coming from the same cell), are combined into one dataset by combining the pts variables of both movies.<br>
<h3>Output</h3><br>
•	For each .HIS file: Matlab data file called _LocRes.mat containing the pts variable of that .HIS file. Example: 20211126_KM12L4a_mEos3.2_03_mov1_LocRes and 20211126_KM12L4a_mEos3.2_03_mov2_LocRes<br>
•	For each .HIS file: PNG file with the plotted reconstructed image with scalebar 5 µm. Example: 20211126_KM12L4a_mEos3.2_03_mov1 and 20211126_KM12L4a_mEos3.2_03_mov2<br>
•	For each .HIS file with the same core name: combined Matlab data file called _Res_all_movies containing the pts variable of those .HIS files. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies<br>
•	For each .HIS file with the same core name: combined PNG file with filename ending in “all” with the plotted reconstructed image with scalebar 5 µm. Example: 20211126_KM12L4a_mEos3.2_03all<br>
<h2>MATLAB_CODE_DBSCAN_and_cluster_properties_bulk<br></h2>
<h3>Input</h3><br>
•	_Res_all_movies file<br>
•	MinVolEllipse function https://nl.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid?requestedDomain= <br>
<h3>Important code aspects<br></h3>
•	Thresholds cell: finds the object with the largest surface area, and only take localizations inside this area into account<br>
•	DBScan analysis with minimum 17 neighbors (MinPts) and search radius (ε) 0.55<br>
•	Calculates cluster properties: the number of the cluster, number of localizations (NrPts), area, density, eccentricity and center of a plotted ellipse (uses MinVolEllipse function)<br>
•	Calculates cell-specific parameters (cell area, number of clusters per cell area, percentage cell area occupied by clusters)<br>
•	Only retains clusters with ≥ 50 localizations and calculates again cell-specific parameters<br>
•	Calculates mean and median cluster specific parameters per cell (both for all clusters and clusters ≥ localizations)<br>
<h3>Output<br></h3>
•	For each _Res_all_movies file: Matlab data file called _Res_all_movies_DBscan containing DBscan variables and results. These include DB, prop_50 and selpts. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies_DBscan<br>
•	For each _Res_all_movies file: Text document with 7 columns of cluster properties, for each of the identified clusters. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies_prop_cl<br>
•	For each _Res_all_movies file: Text document with 7 columns of cluster properties, for each of the identified clusters with ≥ 50 localizations. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies_prop_50<br>
•	For each _Res_all_movies file: Text document with 7 columns of mean and median cluster specific parameters per cell, both for all clusters, and for clusters with  ≥ 50 localizations. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies_means_medians<br>
<h2>MATLAB _CODE_clathrin_cluster_classification_bulk<br></h2>
<h3>Input <br></h3>
•	_Res_all_movies_DBscan<br>
•	MATLAB_CODE_get_prop_clusters<br>
<h3>Important code aspects <br></h3>
•	For every cluster with ≥ 50 localizations, calculate cluster properties. <br>
•	Classify every cluster with ≥ 50 localizations as a pit or a lattice, based on a manually set threshold.<br>
•	Plot the classification result in an image, next to the reconstructed super-resolution image for comparison<br>
•	Per cell, calculate certain cell-specific parameters, including nrFCL per nrPITS<br>
<h3>Output <br></h3>
•	For each _Res_all_movies_DBscan file: Matlab data file called _Res_all_movies_DBscan_th_%th%_individual_clusters_classified. Contains classification results: cellarea; and prop_cl = the individual cluster properties of the classified clusters with ≥ localizations. Example: 20211126_KM12L4a_mEos3.2_03_Res_all_movies_DBscan_th_0.5035_individual_clusters_classified<br>
•	For each _Res_all_movies_DBscan file : a PNG file called _Res_all_movies_DBscan_cluster_classification that shows the classification result next to the super-resolution image. Example: file 20211126_KM12L4a_mEos3.2_03_Res_all_movies_DBscan_cluster_classification<br>
•	Per analyzed file folder a Matlab data file called all_files_cell_classification that contains per analyzed cell the cell-specific parameters in the celldata function.<br>

