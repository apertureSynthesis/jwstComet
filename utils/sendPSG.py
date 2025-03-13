import os,sys

def sendPSG(cfg_file,ret_file,local=True):
    """
    Send a CFG file to the PSG for retrieval and receive the results
    Default choice is using a local copy of PSG instead of the server

    Inputs
        cfg_file - name of the configuration (cfg) file to be sent to the PSG
        ret_file - name of the retrieval/results file where the PSG returned values will be stored
        local - whether we are interrogating a locally installed copy of the PSG or the online version
    """
    if local:
        os.system('curl -d type=ret --data-urlencode file@{} http://localhost:3000/api.php > {}'.format(cfg_file,ret_file))
    else:
        os.system('curl -d type=ret --data-urlencode file@{} https://psg.gsfc.nasa.gov/api.php > {}'.format(cfg_file,ret_file))