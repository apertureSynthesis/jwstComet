import os,sys
import numpy as np
from datetime import datetime
import pandas as pd

"""
Set of functions to facilitate writing and sending configuration (CFG) files for PSG runs
"""

def ephCFG(specFile,name,objectType,midtime,local=True):
    """
    Writes a CFG file that will request the PSG to return a new CFG with ephemeris information. 
    This will be used in the next step to either run a forward model or a retrieval.

    Inputs
        specFile - ASCII file containing the extracted spectrum
        name - name of asteroid or comet
        objectType - type of object, either asteroid or comet
        midtime - midpoint time of the observations
        local - whether we are interrogating a locally installed copy of the PSG

    Outputs
        Relevant configuration files
    """

    #Name of CFG file for requesting the ephemeris parameters
    cfgName = specFile[:-3]+'eph.cfg'
    #Name of CFG file for storing the ephemeris parameters for use in later steps
    ephName = specFile[:-3]+'atm.cfg'

    #Write the CFG file with basic info from the file header
    with open(cfgName, 'w') as fn:
        fn.write('<OBJECT>{}\n'.format(objectType))
        fn.write('<OBJECT-NAME>{}\n'.format(name))
        fn.write('<OBJECT-DATE>{}\n'.format(datetime.strptime(midtime, "%Y-%b-%d %H:%M:%S.%f").strftime("%Y/%m/%d %H:%M")))
        #fn.write('<GEOMETRY>Observatory\n')


    # #Send it to the PSG
    if local:
        os.system('curl -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(cfgName,ephName))
    else:
        os.system('curl -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfgName,ephName))



def atmCFG(specFile, resFile, composition, retrieval, mode, withCont, local=True, masterATM=False, masterATMFile=None):
    """
    Read in the ephemeris CFG file. Incorporate it along with user-defined modeling parameters into a new CFG file.
    Read in the extracted spectrum and header information. Send it off to the PSG for modeling.

    Inputs
        specFile - ASCII file containing the spectrum to be analyzed
        resFile - ASCII file for saving the returned results from the PSG
        composition - dictionary containing compositional information for building the PSG model atmosphere. optional
        retrieval - dictionary containing quantities to be retrieved for each PSG model run. optional.
        mode - extraction mode (circle, rectangle, mapping, azimuthal)
        withCont - whether we are asking the PSG to simulate the continuum or instead simply subtract a baseline

    Outputs
        PSG CFG file containing all parameters, data, and information necessary to run a forward model or retrieval
    """

    #Read in the data file
    wave, spec, err = np.loadtxt(specFile, unpack=1)

    #Dictionary of solar photolysis lifetimes
    solar_lifetimes = {
        'H2O':      {'quiet': 8.294e4, 'active': 4.539e4, 'alias': 'H2O', 'g_alias': 'H2O'},
        'o-H2O':    {'quiet': 8.294e4, 'active': 4.539e4, 'alias': 'H2O_ortho', 'g_alias': 'H2O'},
        'p-H2O':    {'quiet': 8.294e4, 'active': 4.539e4, 'alias': 'H2O_para', 'g_alias': 'H2O'},
        'HDO':      {'quiet': 8.294e4, 'active': 4.539e4, 'alias': 'HDO', 'g_alias': 'H2O'},
        'OHP':      {'quiet': 8.294e4, 'active': 4.539e4, 'alias': 'OHP', 'g_alias': 'OH'},
        'CO2':      {'quiet': 4.948e5, 'active': 2.101e5, 'alias': 'CO2', 'g_alias': 'CO2'},
        '(13)CO2':    {'quiet': 4.948e5, 'active': 2.101e5, 'alias': '13CO2', 'g_alias': 'CO2'},
        'OCS':      {'quiet': 9.807e3, 'active': 7.723e3, 'alias': 'OCS', 'g_alias': 'OCS'},
        'HCN':      {'quiet': 7.662e4, 'active': 3.085e4, 'alias': 'HCN', 'g_alias': 'HCN'},
        'CO':       {'quiet': 1.335e6, 'active': 5.320e5, 'alias': 'CO', 'g_alias': 'CO'},
        '13CO':     {'quiet': 1.335e6, 'active': 5.320e5, 'alias': 'CO_1316', 'g_alias': 'CO'},
        'C17O':     {'quiet': 1.335e6, 'active': 5.320e5, 'alias': 'CO_1217', 'g_alias': 'CO'},
        'C18O':     {'quiet': 1.335e6, 'active': 5.320e5, 'alias': 'CO_1218', 'g_alias': 'CO'},
        'H2CO':     {'quiet': 4.649e3, 'active': 4.369e3, 'alias': 'H2CO', 'g_alias': 'H2CO'},
        'CH4':      {'quiet': 1.317e5, 'active': 5.381e4, 'alias': 'CH4', 'g_alias': 'CH4'},
        'CH3D':     {'quiet': 1.317e5, 'active': 5.381e4, 'alias': 'CH3D', 'g_alias': 'CH4'},
        'C2H6':     {'quiet': 9.491e4, 'active': 3.978e4, 'alias': 'C2H6', 'g_alias': 'C2H6'},
        'CH3OH':    {'quiet': 8.787e4, 'active': 4.816e4, 'alias': 'CH3OH', 'g_alias': 'CH3OH'},
        'CH3OH_V9': {'quiet': 8.787e4, 'active': 4.816e4, 'alias': 'CH3OH_V9', 'g_alias': 'CH3OH'},
        'NH3':      {'quiet': 5.658e3, 'active': 5.022e3, 'alias': 'NH3', 'g_alias': 'NH3'},
        'C2H2':     {'quiet': 3.269e4, 'active': 1.691e4, 'alias': 'C2H2', 'g_alias': 'C2H2'},
        'CN':       {'quiet': '1.3e4 2.1e5', 'active': '1.3e4 2.1e5', 'alias': 'CN', 'g_alias': 'CN'},
        'NH2':      {'quiet': '4.1e3 6.2e4', 'active': '4.1e3 6.2e4', 'alias': 'NH2', 'g_alias': 'NH2'},
        'H2CO-Daughter': {'quiet': '1.428e3 4.649e3', 'active': '1.428e3 4.369e3', 'alias': 'H2CO', 'g_alias': 'H2CO'},
        'NH':       {'quiet': '5.0e4 1.5e5', 'active': '5.0e4 1.5e5', 'alias': 'NH', 'g_alias': 'NH'}
    }

    #Dictionary of spectral resolutions
    resolution = {
        'NIRSPEC':{
            'PRISM/CLEAR':  {'low': 0.60/30., 'high': 5.30/330.},
            'G140M/F070LP': {'low': 0.70/500., 'high': 1.26/898.},
            'G140M/F100LP': {'low': 0.98/699., 'high': 1.88/1343.},
            'G235M/F170LP': {'low': 1.70/722., 'high': 3.15/1342.},
            'G395M/F290LP': {'low': 2.88/728., 'high': 5.20/1317.},
            'G140H/F070LP': {'low': 0.70/1321., 'high': 1.26/2395.},
            'G140H/F100LP': {'low': 0.98/1849., 'high': 1.87/3675.},
            'G235H/F170LP': {'low': 1.70/1911., 'high': 3.15/3690.},
            'G395H/F290LP': {'low': 2.88/1927., 'high': 5.20/3613.}
        },
        'MIRI':{
            '1/SHORT': {'low': 4.885/3400, 'high': 5.751/4000},
            '1/MEDIUM': {'low': 5.634/3420, 'high': 6.632/3990},
            '1/LONG': {'low': 6.408/3330, 'high': 7.524/3840},
            '2/SHORT': {'low': 7.477/3190, 'high': 8.765/3620},
            '2/MEDIUM': {'low': 8.711/3040, 'high': 10.228/3530},
            '2/LONG': {'low': 10.017/2890, 'high': 11.753/3374},
            '3/SHORT': {'low': 11.481/2450, 'high': 13.441/3010},
            '3/MEDIUM': {'low': 13.319/2300, 'high': 15.592/2460},
            '3/LONG': {'low': 15.400/2020, 'high': 18.072/2790},
            '4/SHORT': {'low':  17.651/1400, 'high': 20.93/1960},
            '4/MEDIUM': {'low': 20.417/1660, 'high': 24.220/1730},
            '4/LONG': {'low': 23.884/1340, 'high': 28.329/1520},
            '1/MULTIPLE': {'low': 4.885/3400, 'high': 7.524/3840},
            '2/MULTIPLE': {'low': 7.477/3190, 'high': 11.753/3374},
            '3/MULTIPLE': {'low': 11.481/2450, 'high': 18.072/2790},
            '4/MULTIPLE': {'low': 17.651/1400, 'high': 28.329/1520}
        }
    }


    atm_keys = list(composition.keys())

    gases = atm_keys[9:]
    gas_names = gases.copy()
    for i in range(len(gases)):
            gas_names[i] = solar_lifetimes[gases[i]]['g_alias']

    n_gas = len(gases)

    #Add the retrieval parameters
    ret_keys = list(retrieval.keys())
    if 'RP' in retrieval:
        ret_vars = ret_keys[11:]
    else:
        ret_vars = ret_keys[10:]
    n_vars = len(ret_vars)

    #Read in the CFG file and update it to fit our desired atmosphere and retrieval parameters if any of them are already in the CFG
    #With update as of 4/14/25, the only default parameters are STRUCTURE, PRESSURE, PUNIT, and WEIGHT, so now we have to be sure
    #to add the rest in ourselves.
    retName = specFile[:-3]+'ret.cfg'
    if masterATM:
        os.system('cp {} {}'.format(masterATMFile,retName))
        atmTemplate = masterATMFile
    else:
        os.system('cp {} {}'.format(specFile[:-3]+'atm.cfg',retName))
        atmTemplate = specFile[:-3]+'atm.cfg'

    with open(atmTemplate) as an:
        with open(retName, 'w') as fn:
            for _, line in enumerate(an):
                if '<ATMOSPHERE-PRESSURE>' in line:
                    modified_line = '<ATMOSPHERE-PRESSURE>{}\n'.format(composition['COMA-ACTIVITY']['value'])
                    fn.write(modified_line)
                elif '<ATMOSPHERE-TEMPERATURE>' in line:
                    modified_line = '<ATMOSPHERE-TEMPERATURE>{}\n'.format(composition['TEMPERATURE']['value'])
                    fn.write(modified_line)
                elif '<ATMOSPHERE-NGAS>' in line:
                    modified_line = '<ATMOSPHERE-NGAS>{}\n'.format(n_gas)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-GAS>' in line:
                    gas_list = ','.join(gas_names)
                    modified_line = '<ATMOSPHERE-GAS>{}\n'.format(gas_list)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-TYPE>' in line:
                    models = ['GSFC[' + solar_lifetimes[i]['alias'] + ']' for i in gases]
                    model_list = ','.join(models)
                    modified_line = '<ATMOSPHERE-TYPE>{}\n'.format(model_list)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-ABUN>' in line:
                    abunds = [str(composition[i]['value']) for i in gases]
                    abund_list = ','.join(abunds)
                    modified_line = '<ATMOSPHERE-ABUN>{}\n'.format(abund_list)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-UNIT>' in line:
                    units = [composition[i]['unit'] for i in gases]
                    unit_list = ','.join(units)
                    modified_line = '<ATMOSPHERE-UNIT>{}\n'.format(unit_list)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-TAU>' in line:
                    if composition['Solar Activity'] == 'active':
                        lifetimes = [str(solar_lifetimes[i]['active']) for i in gases]
                    elif composition['Solar Activity'] == 'quiet':
                        lifetimes = [str(solar_lifetimes[i]['quiet']) for i in gases]
                    else:
                        raise ValueError('Must specify active or quiet solar activity levels')
                    
                    lifetime_list = ','.join(lifetimes)
                    modified_line = '<ATMOSPHERE-TAU>{}\n'.format(lifetime_list)
                    fn.write(modified_line)
                elif '<OBJECT-DIAMETER>' in line:
                    continue
                elif '<ATMOSPHERE-WEIGHT>' in line:
                    if composition['ATMOSPHERE-WEIGHT']['set']:
                        vexp = composition['ATMOSPHERE-WEIGHT']['value']
                        modified_line = '<ATMOSPHERE-WEIGHT>{}\n'.format(vexp)
                        fn.write(modified_line)
                    else:
                        fn.write(line)
                else:
                    fn.write(line)
            # #Add in the atmospheric properties
            # modified_line = '<ATMOSPHERE-STRUCTURE>Coma\n'
            # fn.write(modified_line)

            # modified_line = '<ATMOSPHERE-PRESSURE>{}\n'.format(composition['COMA-ACTIVITY']['value'])
            # fn.write(modified_line)

            # modified_line = '<ATMOSPHERE-TEMPERATURE>{}\n'.format(composition['TEMPERATURE']['value'])
            # fn.write(modified_line)

            # modified_line = '<ATMOSPHERE-NGAS>{}\n'.format(n_gas)
            # fn.write(modified_line)

            # gas_list = ','.join(gases)
            # modified_line = '<ATMOSPHERE-GAS>{}\n'.format(gas_list)
            # fn.write(modified_line)

            # models = ['GSFC[' + solar_lifetimes[i]['alias'] + ']' for i in gases]
            # model_list = ','.join(models)
            # modified_line = '<ATMOSPHERE-TYPE>{}\n'.format(model_list)
            # fn.write(modified_line)

            # abunds = [str(composition[i]['value']) for i in gases]
            # abund_list = ','.join(abunds)
            # modified_line = '<ATMOSPHERE-ABUN>{}\n'.format(abund_list)
            # fn.write(modified_line)

            # units = [composition[i]['unit'] for i in gases]
            # unit_list = ','.join(units)
            # modified_line = '<ATMOSPHERE-UNIT>{}\n'.format(unit_list)
            # fn.write(modified_line)

            if composition['Solar Activity'] == 'active':
                lifetimes = [str(solar_lifetimes[i]['active']) for i in gases]
            elif composition['Solar Activity'] == 'quiet':
                lifetimes = [str(solar_lifetimes[i]['quiet']) for i in gases]
            else:
                raise ValueError('Must specify active or quiet solar activity levels')
            
            lifetime_list = ','.join(lifetimes)
            modified_line = '<ATMOSPHERE-TAU>{}\n'.format(lifetime_list)
            fn.write(modified_line)


            #Finish adding the continuum properties
            if retrieval['COMA-OPACITY'] == 'thin':
                fn.write('<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all,FluorThin\n')
            else:
                fn.write('<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all\n')
            modified_line = '<SURFACE-GAS-RATIO>{}\n'.format(composition['SURFACE-GAS-RATIO']['value'])
            fn.write(modified_line)
            modified_line = '<SURFACE-GAS-UNIT>{}\n'.format(composition['SURFACE-GAS-RATIO']['unit'])
            fn.write(modified_line)
            modified_line = '<SURFACE-TEMPERATURE>{}\n'.format(composition['SURFACE-TEMPERATURE']['value'])
            fn.write(modified_line)
            modified_line = '<SURFACE-ALBEDO>{}\n'.format(composition['SURFACE-ALBEDO']['value'])
            fn.write(modified_line)
            modified_line = '<SURFACE-EMISSIVITY>{}\n'.format(composition['SURFACE-EMISSIVITY']['value'])
            fn.write(modified_line)

            #Add in properties for JWST
            fn.write('<GENERATOR-TELESCOPE>SINGLE\n')
            fn.write('<GENERATOR-DIAMTELE>5.64\n')
            fn.write('<GENERATOR-TELESCOPE1>1\n')
            fn.write('<GENERATOR-TELESCOPE2>2.0\n')
            fn.write('<GENERATOR-TELESCOPE3>1.0\n')
            fn.write('<GENERATOR-NOISE>NO\n')

            #Add in flux and modeling properties
            fn.write('<GENERATOR-GAS-MODEL>Y\n')
            if withCont:
                fn.write('<GENERATOR-CONT-MODEL>Y\n')
            else:
                fn.write('<GENERATOR-CONT-MODEL>N\n')
            fn.write('<GENERATOR-LOGRAD>N\n')
            fn.write('<GENERATOR-CONT-STELLAR>Y\n')
            fn.write('<GENERATOR-RADUNITS>Jy\n')
            fn.write('<GENERATOR-TRANS-APPLY>N\n')
            fn.write('<GENERATOR-RESOLUTIONKERNEL>Y\n')
            fn.write('<GENERATOR-TRANS>02-01\n')

            #Add in geometric factors such as offsets and wavelength ranges
            fn.write('<GENERATOR-RANGE1>{}\n'.format(wave[0]))
            fn.write('<GENERATOR-RANGE2>{}\n'.format(wave[-1]))
            with open(specFile) as sf:
                for _, line in enumerate(sf):
                    if '#Instrument used' in line:
                        instrument = line.split()[-1]
                    if '#Grating used' in line:
                        grating = line.split()[-1]
                    if '#Pixel scale (arcsec/pixel)' in line:
                        psa = float(line.split()[-1])
                    if '#Aperture radius (arcsec)' in line:
                        radAp = float(line.split()[-1])
                    if '#X offset (arcsec)' in line:
                        xOffset = float(line.split()[-1])
                    if '#Y offset (arcsec)' in line:
                        yOffset = float(line.split()[-1])
                    if '#Inner annulus radius (arcsec)' in line:
                        innerRadius = float(line.split()[-1])
                    if '#Outer annulus radius (arcsec)' in line:
                        outerRadius = float(line.split()[-1])
                    if '#Aperture size (arcsec)' in line:
                        radWidth = float(line.split()[-2])
                        radHeight = float(line.split()[-1])
            #Work out resolution from dictionary
            if 'RP' in retrieval:
                res_element = retrieval['RP']['res_element']
                res_type = retrieval['RP']['res_type']
            else:
                if instrument == 'NIRSPEC':
                    if grating == 'G395M/F290LP':
                        res_element = np.sqrt(2*np.log(2))*2*0.001440
                        res_type = 'um'
                    elif (grating == 'G395H/F290LP') or (grating == 'G140H/F100LP') or (grating == 'G235H/F170LP'):
                        res_element = 2700
                        res_type = 'RP'
                    elif (grating == 'PRISM/CLEAR'):
                        res_element = 0.022
                        res_type = 'um'
                    # elif grating == 'G395H/F290LP':
                    #     res_element = np.sqrt(2*np.log(2))*2*7.236e-04
                    # elif grating == 'G140H/F100LP':
                    #     res_element = np.sqrt(2*np.log(2))*2*2.388e-04
                    # elif grating == 'G235H/F170LP':
                    #     res_element = np.sqrt(2*np.log(2))*2*3.344e-04
                    else:
                        res_element = 0.5*(resolution[instrument][grating]['low'] + resolution[instrument][grating]['high'])
                        res_type = 'um'
                elif instrument == 'MIRI':
                    if '1/' in grating:
                        res_element = np.sqrt(2*np.log(2))*2*0.000828
                    else:
                        res_element = 0.5*(resolution[instrument][grating]['low'] + resolution[instrument][grating]['high'])
                    res_type = 'um'
            fn.write('<GENERATOR-RESOLUTION>{}\n'.format(res_element))
            fn.write('<GENERATOR-RESOLUTIONUNIT>{}\n'.format(res_type))
            if mode == 'circle':
                fn.write('<GENERATOR-BEAM>{}\n'.format(2*radAp))
                fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-NS>{}\n'.format(yOffset))
                fn.write('<GEOMETRY-OFFSET-EW>{}\n'.format(xOffset))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')
            elif mode == 'rectangle':
                fn.write('<GENERATOR-BEAM>{},{},0,R\n'.format(radWidth,radHeight))
                fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-NS>{}\n'.format(yOffset))
                fn.write('<GEOMETRY-OFFSET-EW>{}\n'.format(xOffset))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')
            elif mode == 'mapping':
                fn.write('<GENERATOR-BEAM>{},{},0,R\n'.format(radWidth,radHeight))
                fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-NS>{}\n'.format(yOffset))
                fn.write('<GEOMETRY-OFFSET-EW>{}\n'.format(xOffset))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')
            elif mode == 'azimuthal':
                if '#Inner annulus radius (arcsec)' in line:
                    innerRadius = float(line.split()[-1])
                if '#Outer annulus radius (arcsec)' in line:
                    outerRadius = float(line.split()[-1])
                if innerRadius == 0:
                    fn.write('<GENERATOR-BEAM>{}\n'.format(2*outerRadius))
                    fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                else:
                    fn.write('<GENERATOR-BEAM>{},{},0,R\n'.format(psa,psa))
                    fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-EW>0\n')
                if innerRadius == 0:
                    fn.write('<GEOMETRY-OFFSET-NS>0\n')
                else:
                    fn.write('<GEOMETRY-OFFSET-NS>{}\n'.format(0.5*(innerRadius + outerRadius)))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')    
            else:
                print('Allowed modes are circle, rectangle, azimuthal, and mapping.')
                return

            fn.write('<OBJECT-DIAMETER>{}\n'.format(composition['DIAMETER']['value']))                

            #Hard code some retrieval preferences, make room to update later
            fn.write('<RETRIEVAL-GAMMA>{}\n'.format(retrieval['GAMMA']))
            fn.write('<RETRIEVAL-ALPHA>\n')
            fn.write('<RETRIEVAL-FLUXSCALER>{}\n'.format(retrieval['FLUXSCALER']))
            fn.write('<RETRIEVAL-FITTELLURIC>{}\n'.format(retrieval['FITTELLURIC']))
            fn.write('<RETRIEVAL-FITGAIN>{}\n'.format(retrieval['FITGAIN']))
            fn.write('<RETRIEVAL-FITGAIN-PHOTOMETRIC>{}\n'.format(retrieval['FITGAIN-PHOTOMETRIC']))
            fn.write('<RETRIEVAL-REMOVEOFFSET>{}\n'.format(retrieval['REMOVEOFFSET']))
            fn.write('<RETRIEVAL-REMOVEFRINGE>{}\n'.format(retrieval['REMOVEFRINGE']))
            fn.write('<RETRIEVAL-FITSTELLAR>{}\n'.format(retrieval['FITSTELLAR']))
            fn.write('<RETRIEVAL-FITFREQ>{}\n'.format(retrieval['FITFREQ']))
            fn.write('<RETRIEVAL-FREQSHIFT>\n')
            fn.write('<RETRIEVAL-NVARS>{}\n'.format(n_vars))
            fn.write('<RETRIEVAL-VARIABLES>{}\n'.format(','.join(ret_vars)))

            fn.write('<RETRIEVAL-VALUES>{}\n'.format(','.join([str(retrieval[i]['start']) for i in ret_vars])))
            fn.write('<RETRIEVAL-MIN>{}\n'.format(','.join([str(retrieval[i]['min']) for i in ret_vars])))
            fn.write('<RETRIEVAL-MAX>{}\n'.format(','.join([str(retrieval[i]['max']) for i in ret_vars])))
            fn.write('<RETRIEVAL-UNITS>{}\n'.format(','.join([retrieval[i]['unit'] for i in ret_vars])))
            fn.write('<DATA>\n')
            for i in range(len(wave)):
                if not np.isnan(spec[i]):
                    fn.write('{} {} {}\n'.format(wave[i],spec[i],err[i]))
                else:
                    continue
            fn.write('</DATA>\n')

    #Send it to the PSG
    if local:
        os.system('curl -d type=ret --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(retName,resFile))
    else:
        os.system('curl -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(retName,resFile))



    

