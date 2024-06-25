import os,sys
import numpy as np
from datetime import datetime
import pandas as pd

def ephCFG(specFile,name,objectType,midtime,delta,key=None):
    """
    Writes a CFG file that will request the PSG to return an updated CFG with ephemeris information. 
    This will be used in the next step to either run a forward model of a retrieval.
    """

    cfgName = specFile[:-3]+'eph.cfg'
    ephName = specFile[:-3]+'atm.cfg'

    #Write the CFG file with basic info
    with open(cfgName, 'w') as fn:
        fn.write('<OBJECT>{}\n'.format(objectType))
        fn.write('<OBJECT-NAME>{}\n'.format(name))
        fn.write('<OBJECT-DATE>{}\n'.format(datetime.strptime(midtime, "%Y-%b-%d %H:%M:%S.%f").strftime("%Y/%m/%d %H:%M")))
        fn.write('<GEOMETRY>Observatory\n')
        #fn.write('<GEOMETRY-OBS-ALTITUDE>{}\n'.format(delta))
        #fn.write('<GEOMETRY-ALTITUDE-UNIT>AU\n')

    #Send it to the PSG
    if key == None:
        #os.system('curl -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfgName,ephName))
        os.system('curl -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(cfgName,ephName))
    else:
        #os.system('curl -d key={} -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(key,cfgName,ephName))
        os.system('curl -d key={} -d type=cfg -d wgeo=y -d wephm=y -d watm=y --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(key,cfgName,ephName))


def atmCFG(specFile, resFile, composition, retrieval, mode, withCont, key=None):

    #Read in the data file
    wave, spec, err = np.loadtxt(specFile, unpack=1)


    #Dictionary of solar photolysis lifetimes
    solar_lifetimes = {
        'H2O':      {'quiet': 8.294e4, 'active': 4.539e4},
        'OHP':      {'quiet': 8.294e4, 'active': 4.539e4},
        'CO2':      {'quiet': 4.948e5, 'active': 2.101e5},
        '13CO2':    {'quiet': 4.948e5, 'active': 2.101e5},
        'OCS':      {'quiet': 9.807e3, 'active': 7.723e3},
        'HCN':      {'quiet': 7.662e4, 'active': 3.085e4},
        'CO':       {'quiet': 1.335e6, 'active': 5.320e5},
        'H2CO':     {'quiet': 4.649e3, 'active': 4.369e3},
        'CH4':      {'quiet': 1.317e5, 'active': 5.381e4},
        'C2H6':     {'quiet': 9.491e4, 'active': 3.978e4},
        'CH3OH':    {'quiet': 8.787e4, 'active': 4.816e4},
        'CH3OH_V9': {'quiet': 8.787e4, 'active': 4.816e4},
        'NH3':      {'quiet': 5.658e3, 'active': 5.022e3},
        'C2H2':     {'quiet': 3.269e4, 'active': 1.691e4},
        'CN':       {'quiet': '1.3e4 2.1e5', 'active': '1.3e4 2.1e5'},
        'NH2':      {'quiet': '4.1e3 6.2e4', 'active': '4.1e3 6.2e4'}
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
            '1/SHORT': 0.0009,
            '1/MEDIUM': 0.0011,
            '1/LONG': 0.0012,
            '2/SHORT': 0.0014,
            '2/MEDIUM': 0.0016,
            '2/LONG': 0.0019,
            '3/SHORT': 0.0021,
            '3/MEDIUM': 0.0025,
            '3/LONG': 0.0029,
            '4/SHORT': 0.0035,
            '4/MEDIUM': 0.0042,
            '4/LONG': 0.0049,
            '1/MULTIPLE': 0.00105,
            '2/MULTIPLE': 0.00165,
            '3/MULTIPLE': 0.0025,
            '4/MULTIPLE': 0.0042
        }
    }


    atm_keys = list(composition.keys())

    gases = atm_keys[8:]
    n_gas = len(gases)

    #Add the retrieval parameters
    ret_keys = list(retrieval.keys())
    ret_vars = ret_keys[10:]
    n_vars = len(ret_vars)

    #Read in the CFG file and update it to fit our desired atmosphere and retrieval parameters
    retName = specFile[:-3]+'ret.cfg'
    os.system('cp {} {}'.format(specFile[:-3]+'atm.cfg',retName))
    with open(specFile[:-3]+'atm.cfg') as an:
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
                    gas_list = ','.join(gases)
                    modified_line = '<ATMOSPHERE-GAS>{}\n'.format(gas_list)
                    fn.write(modified_line)
                elif '<ATMOSPHERE-TYPE>' in line:
                    models = ['GSFC[' + i + ']' for i in gases]
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
                else:
                    fn.write(line)

            #Finish adding the continuum properties
            if withCont:
                fn.write('<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all\n')
                modified_line = '<SURFACE-GAS-RATIO>{}\n'.format(composition['SURFACE-GAS-RATIO']['value'])
                fn.write(modified_line)
                modified_line = '<SURFACE-GAS-UNIT>{}\n'.format(composition['SURFACE-GAS-RATIO']['value'])
                fn.write(modified_line)
                modified_line = '<SURFACE-TEMPERATURE>{}\n'.format(composition['SURFACE-TEMPERATURE']['value'])
                fn.write(modified_line)
                modified_line = '<SURFACE-ALBEDO>{}\n'.format(composition['SURFACE-ALBEDO']['value'])
                fn.write(modified_line)
                modified_line = '<SURFACE-EMISSIVITY>{}\n'.format(composition['SURFACE-EMISSIVITY']['value'])
                fn.write(modified_line)
                fn.write('<GENERATOR-CONT-MODEL>Y\n')
            else:
                fn.write('<GENERATOR-CONT-MODEL>N\n')

            #Add in properties for JWST
            fn.write('<GENERATOR-TELESCOPE>SINGLE\n')
            fn.write('<GENERATOR-DIAMTELE>5.64\n')

            #Add in flux and modeling properties
            fn.write('<GENERATOR-GAS-MODEL>Y\n')
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
                    # if '#Spectral plate scale (um)' in line:
                    #     res_element = float(line.split()[-1])
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
            #Work out resolution from dictionary
            if instrument == 'NIRSPEC':
                res_element = 0.5*(resolution[instrument][grating]['low'] + resolution[instrument][grating]['high'])
            elif instrument == 'MIRI':
                res_element = resolution[instrument][grating]
            fn.write('<GENERATOR-RESOLUTION>{}\n'.format(res_element))
            fn.write('<GENERATOR-RESOLUTIONUNIT>um\n')
            if mode == 'beam':
                fn.write('<GENERATOR-BEAM>{}\n'.format(2*radAp))
                fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-NS>{}\n'.format(yOffset))
                fn.write('<GEOMETRY-OFFSET-EW>{}\n'.format(xOffset))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')
            if mode == 'azimuthal':
                if '#Pixel scale (arcsec/pixel)' in line:
                    psa = float(line.split()[-1])
                if '#Inner annulus radius (arcsec)' in line:
                    innerRadius = float(line.split()[-1])
                if '#Outer annulus radius (arcsec)' in line:
                    outerRadius = float(line.split()[-1])
                if innerRadius == 0:
                    fn.write('<GENERATOR-BEAM>{}\n'.format(2*outerRadius))
                    fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                else:
                    fn.write('<GENERATOR-BEAM>{},{},0,R\n'.format(psa,outerRadius - innerRadius))
                    fn.write('<GENERATOR-BEAM-UNIT>arcsec\n')
                fn.write('<GEOMETRY-OFFSET-NS>0\n')
                if innerRadius == 0:
                    fn.write('<GEOMETRY-OFFSET-EW>0\n')
                else:
                    fn.write('<GEOMETRY-OFFSET-EW>{}\n'.format(0.5*(innerRadius + outerRadius)))
                fn.write('<GEOMETRY-OFFSET-UNIT>arcsec\n')    

            fn.write('<OBJECT-DIAMETER>{}\n'.format(composition['DIAMETER']['value']))                

            #Hard code some retrieval preferences, make room to update later
            fn.write('<RETRIEVAL-GAMMA>{}\n'.format(retrieval['GAMMA']))
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
                    if retrieval['COMA-OPACITY'] == 'thin':
                        fn.write('{} {} {}\n'.format(wave[i],spec[i]/1000.,err[i]/1000.))
                    else:
                        fn.write('{} {} {}\n'.format(wave[i],spec[i],err[i]))
                else:
                    continue
            fn.write('</DATA>\n')

    #Send it to the PSG
    if key == None:
        #os.system('curl -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(retName,resFile))
        os.system('curl -d type=ret --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(retName,resFile))
    else:
        #os.system('curl -d key={} -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(key,retName,resFile))
        os.system('curl -d key={} -d type=ret --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(key,retName,resFile))



    

