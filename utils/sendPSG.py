import os,sys

def sendPSG(cfg_file,ret_file,key=None):
    """
    Send a CFG file to the PSG for retrieval and receive the results
    """
    if key != None:
        os.system('curl -d key={} -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(key,cfg_file,ret_file))
    else:
        os.system('curl -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfg_file,ret_file))