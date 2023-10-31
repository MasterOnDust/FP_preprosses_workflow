
SDATE=config['sdate']
EDATE=config['edate']

SDATE_m = config['m_sdate']
EDATE_m = config['m_edate']

rule calc_receptor_correlations:
    input:
        emission_data = expand(config['flexdust_results']+'/emission_flux.time_series.{region}.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                        region=['taklamakan','mongolia','north_west', 'total'],allow_missing=True),
        receptor_data_wetdep_2micron = expand('results/model_results/time_series/wetdep/wetdep.{location}.{source}.2micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys() , 
                               source=['mongolia','taklamakan','north_west', 'total']),
        receptor_data_drydep_2micron = expand('results/model_results/time_series/drydep/drydep.{location}.{source}.2micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys(),
                               source=['mongolia','taklamakan','north_west', 'total']),
        receptor_data_wetdep_20micron = expand('results/model_results/time_series/wetdep/wetdep.{location}.{source}.20micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys(),
                               source=['mongolia','taklamakan','north_west', 'total']),
        receptor_data_drydep_20micron = expand('results/model_results/time_series/drydep/drydep.{location}.{source}.20micron.MAM.'+ str(SDATE_m)+'-'+str(EDATE_m)+'.nc',
                               location=config['receptors'].keys(),
                               source=['mongolia','taklamakan','north_west', 'total'])
        
    output:
        outpath='results/timeseries_table.csv'
    run:
        import pandas as pd
        import xarray as xr
        from thesis_toolbox.composites.create_composites import detrend_timeseries
        def read_deposition_data(paths, psize, detrend=False):
            dsets = []
            for path in paths:
                receptor_name = path.split('/')[-1].split('.')[1]
                kind = path.split('/')[-1].split('.')[0]
                source_reg = path.split('/')[-1].split('.')[2]
                ds = xr.open_dataset(path)
                vname = f'{source_reg} {receptor_name} {kind} {psize}'
                ds = ds.rename({ds.varName:vname})
                if detrend:
                    dsets.append(detrend_timeseries(ds[vname].to_series()))
                else:
                    dsets.append(ds[vname].to_series())
            
            return dsets

        taklamakan =  xr.open_dataset(input.emission_data[0])
        mongolia = xr.open_dataset(input.emission_data[1])
        north_west = xr.open_dataset(input.emission_data[2])
        total = xr.open_dataset(input.emission_data[3])
        
        emission_data = {
            'Emissions taklamakan' : taklamakan[taklamakan.varName].values,
            'Emissions mongolia' : mongolia[mongolia.varName].values,
            'Emissions north_west' : north_west[north_west.varName].values,
            'Emissions total' : total[total.varName].values
        }

        df_emission_data = pd.DataFrame(emission_data, index=taklamakan.time.dt.year.values)
        data = df_emission_data
        deposition_data = {
            'wet_2micron' : read_deposition_data(input.receptor_data_wetdep_2micron,'2micron', detrend=False),
            'wet_20micron' : read_deposition_data(input.receptor_data_wetdep_20micron,'20micron',detrend=False),
            'dry_2micron' :  read_deposition_data(input.receptor_data_drydep_2micron,'2micron',detrend=False),
            'dry_20micron' : read_deposition_data(input.receptor_data_drydep_20micron,'20micron',detrend=False)  
        }
        
        dfs = []
        for key,depodata in deposition_data.items():
            dfs.append(
                pd.DataFrame(depodata).T
            )
        data = data.join(dfs)
        data.to_csv(output.outpath)

