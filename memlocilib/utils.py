from . import config



def get_biopy_pssm(sequence, profile_matrix):
    from Bio.Align import AlignInfo
    alph = "ARNDCQEGHILKMFPSTWYV"
    biopy_pssm = []
    for i in range(len(sequence)):
        biopy_pssm.append((sequence[i], {alph[j]:profile[i][j] for j in range(len(alph))}))
    return AlignInfo.PSSM(biopy_pssm)

def cut_peptide(i_json):
    peptide = [f for f in i_json['features'] if f['type'] == "SIGNAL" or f['type'] == "TRANSIT"]
    cleavage = 0
    sequence = i_json['sequence']['sequence']
    if len(peptide) > 0:
        cleavage = peptide[0]['end']
        sequence = i_json['sequence']['sequence'][cleavage:]
    return sequence, cleavage

def get_json_output(i_json, memloci_pred):
    loc = memloci_pred[1]
    score = float(memloci_pred[2][loc][:-1])/100.0
    if 'comments' not in i_json:
        i_json['comments'] = []
    if 'dbReferences' not in i_json:
        i_json['dbReferences'] = []

    go_info = config.GOINFO[loc]
    i_json['dbReferences'].append({
        "id": go_info['goid'],
        "type": "GO",
        "properties": {
          "term": go_info['term'],
          "source": "IEA:MemLoci",
          "score": round(float(score),2)
        },
        "evidences": [
          {
            "code": "ECO:0000256",
            "source": {
              "name": "SAM",
              "id": "MemLoci",
              "url": "https://mu2py.biocomp.unibo.it/memloci/",
            }
          }
        ]
    })
    if len(i_json['comments']) == 0:
        i_json['comments'].append({
            "type": "SUBCELLULAR_LOCATION",
            "locations": [
              {
                "location": {
                  "value": go_info["uniprot"],
                  "score": round(float(score),2),
                  "evidences": [
                    {
                      "code": "ECO:0000256",
                      "source": {
                        "name": "SAM",
                        "id": "MemLoci",
                        "url": "https://mu2py.biocomp.unibo.it/memloci/",
                      }
                    }
                  ]
                }
              }
            ]
        })
    else:
        sl = [c for s in i_json['comments'] if c['type'] == "SUBCELLULAR_LOCATION"][0]
        sl['locations'].append({
          "location": {
            "value": go_info["uniprot"],
            "score": round(float(score),2),
            "evidences": [
              {
                "code": "ECO:0000256",
                "source": {
                  "name": "SAM",
                  "id": "MemLoci",
                  "url": "https://mu2py.biocomp.unibo.it/memloci/",
                }
              }
            ]
          }
        })
        i_json['comments'] = [sl]
    return i_json
