import re
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
import argparse


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--genenames", help="input file containing official gene names from HGNC")
argParser.add_argument("-o", "--output", help="output file")

args = argParser.parse_args()

'''
Output filename
'''
output_filename = args.output

'''
Load HGNC data and explode the Uniprot multiple accessions
'''
hgnc_filename = args.genenames
hgnc_df = pd.read_csv(hgnc_filename, sep='\t')
hgnc_df.rename(columns= {
    'Approved symbol': 'geneSymbol', 
    'UniProt ID(supplied by UniProt)': 'primaryAccession', 
    'Approved name' : 'geneName',
    'HGNC ID': 'HgncId', 
    'Alias symbols': 'AliasSymbols',
    'Ensembl gene ID': 'EnsemblGeneId'}, inplace=True)

# Using Series.astype() to convert to string 
hgnc_df["primaryAccession"]=hgnc_df["primaryAccession"].values.astype(str)
# There may be multiple Uniprot accession separated by comma.
hgnc_df['primaryAccession'] = hgnc_df['primaryAccession'].map(lambda a: a.split(', '))
# convert back to scalar 
hgnc_df["primaryAccession"]=hgnc_df["primaryAccession"].values.astype('object')
hgnc_df = hgnc_df.explode(['primaryAccession'])

'''
Prepare Uniprot API access session
'''
re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25,
                status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        json = response.json()
        total = response.headers["x-total-results"]
        yield json, total
        batch_url = get_next_link(response.headers)


# url = 'https://rest.uniprot.org/uniprotkb/search?fields=accession%2Ccc_interaction&format=tsv&query=Insulin%20AND%20%28reviewed%3Atrue%29&size=500'
url = 'https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606&fields=accession,id,gene_names,ft_intramem,cc_subcellular_location,ft_topo_dom&size=500'
count = 0
total_count = 0
hr = None

# Agree with CI team to use tqdm in the future in place of custom function
printProgressBar(0, 100, prefix = 'Progress:', suffix = 'Complete', length = 20)

for batch, total in get_batch(url):

    data_dict = {
        'primaryAccession': [],
        'uniProtkbId': [],
        'HgncId' : [],
        'geneSymbol': [],
        'synonyms': [],
        'location': [],
        'topology': []
    }

    # Create the main dataframe if it does not exist
    if hr is None:
        hr = pd.DataFrame(data_dict)

    # Loop over the Uniprot records
    for data in batch['results']:
        count += 1
        #print(data["primaryAccession"], data["uniProtkbId"])
        data_dict["primaryAccession"].append(data["primaryAccession"])
        data_dict["uniProtkbId"].append(data["uniProtkbId"])
        geneSymbol = ''
        hgncId = ''
        synonyms = []
        if "genes" in data.keys():
            for gene in data["genes"]:
                if "geneName" in gene.keys():
                    geneSymbol = gene["geneName"]["value"]
                    if 'synonyms' in gene.keys():
                        synonyms = list(
                            map(lambda x: x["value"], gene["synonyms"]))
                    if 'evidences' in gene["geneName"]:
                        for evidence in gene["geneName"]["evidences"]:
                            if "source" in evidence and evidence["source"] == "HGNC":
                                hgncId = evidence["id"]
                                break
                    break
        data_dict["geneSymbol"].append(geneSymbol)
        data_dict["HgncId"].append(hgncId)
        data_dict["synonyms"].append(','.join(synonyms))

        d = {
            'location': [],
            'topology': []
        }

        if "comments" in data:
            for comment in data["comments"]:
                if comment["commentType"] == "SUBCELLULAR LOCATION" and "subcellularLocations" in comment.keys():
                    for subcellularLocation in comment["subcellularLocations"]:
                        for type in ['location', 'topology']:
                            if type in subcellularLocation.keys():
                                value = subcellularLocation[type]["value"]
                                d[type].append(value)
                                # if "evidences" in subcellularLocation[type]:
                                #    for evidence in subcellularLocation[type]["evidences"]:
                                #        value = evidence["evidenceCode"]
                                #        print('evidenceCode', value)
        for type in ['location', 'topology']:
            data_dict[type].append(','.join(set(d[type])))

    # convert to DataFrame (all arrays should be of the same length)
    hr1 = pd.DataFrame(data_dict)

    # add the new rows to the main dataframe
    hr = pd.concat([hr, hr1], ignore_index=True)

    printProgressBar(count, int(total), prefix = 'Progress:', suffix = 'Complete', length = 20)
    
'''
Join HGNC information to add Ensembl identifiers. The majority of HGNC will match
the Uniprot accession.
'''
hgnc_df1 = hgnc_df.drop(['HgncId', 'geneSymbol'], axis=1)
result_uniprot_accession = pd.merge(hr, hgnc_df1, how="left", on=["primaryAccession"])

result_uniprot_accession_filtered = result_uniprot_accession[pd.notnull(result_uniprot_accession['geneName']) | pd.notnull(result_uniprot_accession['EnsemblGeneId'])]

'''
There will be rows with no geneName and no Ensembl id so we subset those and 
join on the gene symbol
'''
missing_df = result_uniprot_accession.loc[pd.isnull(result_uniprot_accession['geneName']) & pd.isnull(result_uniprot_accession['EnsemblGeneId'])]
missing_df = missing_df[['primaryAccession', 'uniProtkbId', 'HgncId', 'geneSymbol', 'synonyms', 'location', 'topology']]
hgnc_df2 = hgnc_df.drop(['primaryAccession', 'HgncId'], axis=1)
result_gene_symbol = pd.merge(missing_df, hgnc_df2, how="left", on=["geneSymbol"])

'''
Now we can concatenate the 2 dataframes
'''
results = pd.concat([result_uniprot_accession_filtered, result_gene_symbol], ignore_index=True, sort=False)

'''
And finally save the file
'''
results.to_csv(output_filename, sep='\t', na_rep='', header=True, index=False)

# Example entry from UniProt API
entry = {"primaryAccession": "P04439",
         "uniProtkbId": "HLAA_HUMAN",
         "genes": [
             {
                 "geneName":
                 {
                     "evidences":
                     [
                        {
                            "evidenceCode": "ECO:0000312",
                            "source": "HGNC",
                            "id": "HGNC:4931"
                        }
                     ],
                     "value": "HLA-A"
                 },
                 "synonyms":
                 [
                     {
                         "value": "HLAA"
                     }
                 ]
             }
         ],
         "comments":
         [
             {
                 "commentType": "SUBCELLULAR LOCATION",
                 "subcellularLocations":
                 [
                     {
                         "location":
                         {
                             "evidences":
                             [
                                {
                                    "evidenceCode": "ECO:0000269",
                                    "source": "PubMed",
                                    "id": "21263072"
                                },
                                 {
                                    "evidenceCode": "ECO:0000269",
                                    "source": "PubMed",
                                    "id": "25880248"
                                },
                                 {
                                    "evidenceCode": "ECO:0000269",
                                    "source": "PubMed",
                                    "id": "8805302"
                                }
                             ],
                             "value": "Cell membrane",
                             "id": "SL-0039"
                         },
                         "topology":
                         {
                             "evidences":
                             [
                                 {
                                     "evidenceCode": "ECO:0000255"
                                 }
                             ],
                             "value": "Single-pass type I membrane protein",
                             "id": "SL-9905"
                         }
                     },
                     {
                         "location":
                         {
                             "evidences": [
                                 {
                                     "evidenceCode": "ECO:0000305",
                                     "source": "PubMed",
                                     "id": "8805302"
                                 }
                             ],
                             "value": "Endoplasmic reticulum membrane",
                             "id": "SL-0097"
                         },
                             "topology":
                             {
                             "evidences":
                             [
                                 {
                                     "evidenceCode": "ECO:0000255"
                                 }
                             ],
                             "value": "Single-pass type I membrane protein",
                             "id": "SL-9905"
                         }
                     }
                 ]
             }
         ],
         "features":
         [
             {
                 "type": "Topological domain",
                 "location":
                 {
                     "start":
                     {
                         "value": 25,
                         "modifier": "EXACT"
                     },
                     "end":
                     {
                         "value": 308,
                         "modifier": "EXACT"
                     }
                 },
                 "description": "Extracellular",
                 "evidences":
                 [
                     {"evidenceCode": "ECO:0000255"
                      }
                 ]
             },
             {
                 "type": "Topological domain",
                 "location":
                 {
                     "start":
                     {
                         "value": 333,
                         "modifier": "EXACT"
                     },
                     "end":
                     {
                         "value": 365, "modifier": "EXACT"
                     }
                 },
                 "description": "Cytoplasmic",
                 "evidences":
                 [
                     {
                         "evidenceCode": "ECO:0000255"
                     }
                 ]
             }
         ]
         }
