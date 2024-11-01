import os,sys

def sendPSG(cfg_file,ret_file,local=True):
    """
    Send a CFG file to the PSG for retrieval and receive the results
    Default choice is using a local copy of PSG instead of the server
    """
    if local:
        os.system('curl -d type=ret --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(cfg_file,ret_file))
    else:
        os.system('curl -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfg_file,ret_file))