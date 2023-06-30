import json
import re
from typing import Dict

# universal functions


def parseGFFAttributes(attributes: str):
    """
    Parse the attribute column of gff file
    """
    parsed = {}
    temp = attributes.split(";")
    for item in temp:
        parsed[item.split("=")[0]] = item.split("=")[1].replace("\n", "")
    return parsed


def jsonString(dict: dict):
    """
    Beautify json
    """
    return json.dumps(dict, indent=4, separators=(",", ": "))


def parseNCBIFastaDescription(header: str) -> Dict[str, str]:
    """
    Convert NCBI Fasta Header to dict

    input example:
        ">accession [key1=value] [key2=value] [key3=value] description"
    """

    data = {}
    for match in re.finditer(r"\[(.*?)\]", header):
        try:
            assert match is not None
            matched_seq = match.group(1)
            temp = matched_seq.split("=")
            assert len(temp) == 2
            key = temp[0]
            assert key not in data.keys()
            value = temp[1]
            data[key] = value
        except BaseException:
            raise ValueError("Description Invalid")
    return data

    """# not in use
def reformatNewickTree(tree: str) -> str:
    
    # Convert multi-branch newick tree to binary tree
    
    results = tree
    for match in re.finditer("\(([^\(\)]+)\)", tree):
        leafs = match.group(1)
        replace = ""
        for index, item in enumerate(leafs.split(",")):
            if index == 0:
                replace = item
            elif index == 1:
                replace = f"{replace},{item}"
            else:
                replace = f"({replace}):0.0,{item}"
        results = results.replace(leafs, replace)
    return results

# Not in use
def convertpdSeries2Dict(row):
    ""
    Convert pandas series to dict
    ""
    info = row.to_dict()
    for key, value in info.items():
        if pd.isna(value) or str(value) == "NA" or str(value) == "NaN" or str(value) == "nan":
            info[key] = None
    return info

    """
