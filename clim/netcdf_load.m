 function [data] = netcdf_load(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads all variables from a netcdf file into a structure
% Usage: data = netcdf_load('dir','./','file','filename.nc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A.dir  = './';
 A.file = 'none';
 A.nan = nan;
% Parse required variables, substituting defaults where necessary
A = parse_pv_pairs(A, varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 tf = fliplr(A.file);
 if ~strcmp(tf(1:3),'cn.');
    A.file = [A.file '.nc'];
 end

 if strcmp(A.file,'none')
    error([A.file ' not found']);
 end

 indir = A.dir;
 fname = A.file;

 try
    ncid = netcdf.open([indir fname],'NC_NOWRITE');
 catch
    error([indir fname '    not found']);
 end

 [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

 for indv=1:numvars
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,indv-1);
    varid = netcdf.inqVarID(ncid,varname);
    varname = lower(varname);
    thisvar = double(netcdf.getVar(ncid,varid)); 
    thisvar(thisvar<=-1e10) = nan;
    data.(varname) = thisvar; 
    if ~isnan(A.nan)
       data.(varname)(data.(varname)==A.nan) = nan;
    end
 end

 netcdf.close(ncid)
 
