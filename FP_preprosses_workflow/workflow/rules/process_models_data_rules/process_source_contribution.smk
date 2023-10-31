

def get_paths(path0,sdate, edate, kind, location, size,ext):
    years = [str(year) for year in range(int(sdate), int(edate)+1)]
    paths=[]
    if size=='2micron':
        tag='CLAY'
    else:
        tag='SILT'
    for year in years:
        temp_path=path0+'/'+kind+'/'+size+'/'+year+'/'+kind
        for dstart,dend in zip(['0301','0331','0430'],['0331', '0430','0531']):
            paths.append(temp_path+'_'+location+'_'+tag+'_'+year+dstart+'-'+year+dend+'.'+ext)

    return paths

rule resample_source_contrib:
    input:
        paths= lambda wildcards: get_paths(config['source_contrib_path'], wildcards.sdate, wildcards.edate,
                                           wildcards.kind, wildcards.location, wildcards.size,'zarr')
        
    output:
        outpath=config['intermediate_results_models']+'/{kind}/{kind}.{location}.{size}.monthly.{sdate}-{edate}.nc'
    wildcard_constraints:
        location='|'.join(config['receptors'].keys()),
        size='2micron|20micron',
        kind='drydep|wetdep',

    threads: 1
    run: 
        from dust.utils.resample import resample_monthly, concatenate_monthly
        import time
        paths=input.paths
        if paths[0].endswith('.zarr'):
            dsets_list=[resample_monthly(xr.open_zarr(path)) for path in paths]
        else:
            dsets_list=[resample_monthly(xr.open_dataset(path)) for path in paths]
        dsets=concatenate_monthly(dsets_list)
        dsets.attrs['history'] = '{} {} '.format(time.ctime(time.time()),'resample_source_contrib') + dsets.attrs['history']
        dsets.to_netcdf(output.outpath)

def get_flexpart_input_paths(w, edate):
    first_piece = f'/{w.kind}/{w.size}/{w.year}/{w.year}' 
    sec_piece = f'_00{w.location}/output'
    path_folder = expand(config['flexpart_path']+first_piece+'{edate}'+sec_piece, edate=edate)
    paths = []
    for folder in path_folder:
        ncfile = glob_wildcards(folder+'/{ncfile}.nc')
        paths.append(folder+'/{}.nc'.format(''.join(ncfile.ncfile))) 
    return paths

rule process_flexpart_output_march:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0331')
    output: 
        outpath=config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0301-{year}0331.nc',
    wildcard_constraints:
        tag='CLAY|SILT'
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1
    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"

rule process_flexpart_output_april:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0430')
    output: 
        outpath=config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0331-{year}0430.nc'
    wildcard_constraints:
        tag='CLAY|SILT'
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1
    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"



rule process_flexpart_output_may:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0531')
    output: 
        config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0430-{year}0531.nc',
    wildcard_constraints:
        tag='CLAY|SILT'
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1

    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"


rule process_flexpart_output_march_zarr:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0331')
    output: 
        outpath=directory(config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0301-{year}0331.zarr')
    wildcard_constraints:
        tag='CLAY|SILT'
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1
    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"

rule process_flexpart_output_april_zarr:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0430')
    output: 
        outpath=directory(config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0331-{year}0430.zarr')
    wildcard_constraints:
        tag='CLAY|SILT',
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1
    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"



rule process_flexpart_output_may_zarr:
    input:
        flexdust_path = config['flexdust_path']+'/{year}/',
        flexpart_path = lambda wildcards: get_flexpart_input_paths(wildcards,edate='0531')
    output: 
        outpath=directory(config['source_contrib_path']+'/{kind}/{size}/{year}/{kind}_{location}_{tag}_{year}0430-{year}0531.zarr')
    wildcard_constraints:
        tag='CLAY|SILT',
    resources:
        time="00:30:00",
        memory_per_job="16GB",
        max_threads=4
    threads:1

    params:
        x0 = config['domain']['lon0'],
        x1 = config['domain']['lon1'],
        y0 = config['domain']['lat0'],
        y1 = config['domain']['lat1'],
        rp = config['source_contrib_path'],
        use_slurm = True,
        use_dask = True
    
    notebook:
        "../../notebooks/process_raw_output.py.ipynb"


rule create_deposition_timeseries:
    input:
        paths= lambda wildcards: get_paths(config['source_contrib_path'], wildcards.year, wildcards.year,
                                           wildcards.kind, wildcards.location, wildcards.size, 'zarr')
        
    output:
        outpath=config['intermediate_results_models']+'/timeseries/{kind}/{kind}.{location}.{size}.{freq}.{year}.csv'
    wildcard_constraints:
        location='|'.join(config['receptors'].keys()),
        size='2micron|20micron',
        kind='drydep|wetdep',
        freq='3H|Day|Weekly'
    run:
        import xarray as xr
        paths = input.paths
        if wildcards.freq=='Day':
            freq='D'
        elif wildcards.freq=='Weekly':
            freq='W'
        else:
            freq=wildcards.freq
        def pre_process(dset):
            dset = dset.drop_vars('surface_sensitivity')
            # dset = dset.drop_dims('numpoint')
            with xr.set_options(keep_attrs=True):
                dset = dset.sum(dim=['btime','lon','lat'])
            return dset
        
        ds = xr.open_mfdataset(paths,concat_dim=['time'], parallel=True, 
                chunks={'time':40}, combine='nested',preprocess=pre_process, engine='zarr')
        if ds.ind_receptor==4 or ds.ind_receptor==3:
            da_dep = ds[wildcards.kind].resample(time=freq).sum()
        else:
            da_dep = ds[wildcards.kind].resample(time=freq).mean()
        df = da_dep.to_pandas()
        df.to_csv(output.outpath)