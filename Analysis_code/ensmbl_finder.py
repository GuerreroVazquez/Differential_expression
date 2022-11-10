import json

import requests, sys



server = "http://rest.ensembl.org"

def get_data_from_ensmbl_id(ensmbl_id):
    ext = f"/lookup/id/{ensmbl_id}?expand=0"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    return repr(decoded)
def get_gene_name_from_ensmbl(gene_id="ENSG00000108861"):
    try:
        x = get_data_from_ensmbl_id(gene_id).replace("\'", "\"")
        x = json.loads(x)
        display_name =  x.get("display_name")
        if display_name is None:
            return 'NF'
    except:
        display_name="E"
    return display_name