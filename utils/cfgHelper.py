import os,sys
from datetime import datetime

def ephCFG(specFile,name,objectType,midtime,delta):
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
        fn.write('<GEOMETRY-OBS-ALTITUDE>{}\n'.format(delta))
        fn.write('<GEOMETRY-ALTITUDE-UNIT>AU\n')

    #Send it to the PSG
    os.system('curl -d type=cfg -d wgeo=y -d wephm=y -d watm=y -d wcon=y --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfgName,ephName))

def atmCFG(specFile, composition, retrieval):

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
        'CN':       {'quiet': '1.3e4 2.1e5', 'active': '1.3e4 2.1e5'}
    }



    atm_keys = list(composition.keys())

    gases = atm_keys[2:]
    n_gas = len(gases)

    #Read in the CFG file and update it to fit our desired atmosphere and retrieval parameters
    retName = specFile[:-3]+'ret.cfg'
    os.system('cp {} {}'.format(specFile[:-3]+'atm.cfg',retName))
    with open(retName, 'w') as fn:
        for _, line in enumerate(fn):
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
                fn.write_modified_line
            elif '<ATMOSPHERE-ABUN>' in line:
                abunds = [composition[i]['value'] for i in gases]
                abund_list = ','.join(abunds)
                modified_line = '<ATMOSPHERE-ABUN>{}\n'.format(abund_list)
                fn.write(modified_line)
            elif '<ATMOSPHERE-UNIT>' in line:
                units = [composition[i]['unit'] for i in gases]
                unit_list = ','.join(units)
                modified_line = '<ATMOSPHERE-UNIT>{}\n'.format(unit_list)
            elif '<ATMOSPHERE-TAU>' in line:
                if composition['Solar activity'] == 'active':
                    lifetimes = [solar_lifetimes[i]['active'] for i in gases]
                elif composition['Solar activity'] == 'quiet':
                    lifetimes = [solar_lifetimes[i]['quiet'] for i in gases]
                else:
                    raise ValueError('Must specify active or quiet solar activity levels')
                
                lifetime_list = ','.join(lifetimes)
                modified_line = '<ATMOSPHERE-TAU>{}\n'.format(lifetime_list)
                fn.write(modified_line)
            else:
                fn.write(line)

            #Finish adding the continuum properties
            fn.write('<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all\n')
            fn.write('<SURFACE-GAS-UNIT>ratio\n')

            #Add the retrieval parameters
            


    

