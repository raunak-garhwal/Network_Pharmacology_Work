import requests

def get_smiles_from_pubchem(compound_name: str):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    # Try name-based property search
    url = f"{base_url}/name/{compound_name}/property/CanonicalSMILES/json"
    r = requests.get(url, timeout=10)

    if r.status_code == 200:
        data = r.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props and "CanonicalSMILES" in props[0]:
            return props[0]["CanonicalSMILES"]

    # Fallback: get CID first
    url = f"{base_url}/name/{compound_name}/cids/json"
    r = requests.get(url, timeout=10)
    if r.status_code == 200:
        data = r.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            cid = cids[0]
            url = f"{base_url}/cid/{cid}/property/CanonicalSMILES/json"
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                props = data.get("PropertyTable", {}).get("Properties", [])
                if props and "CanonicalSMILES" in props[0]:
                    return props[0]["CanonicalSMILES"]

    return None


# Example usage
print(get_smiles_from_pubchem("beta-Bisabolene"))
# print(get_smiles_from_pubchem("eugenol"))