import os
import shutil
from ncbi.datasets.openapi.exceptions import ApiAttributeError
from ncbi.datasets import GenomeApi
from ncbi.datasets.openapi.model.v1_assembly_dataset_request import V1AssemblyDatasetRequest
from ncbi.datasets.openapi.model.v1_annotation_for_assembly_type import V1AnnotationForAssemblyType
import subprocess

"""
    Encapsulation of ncbi apis for downloading assemblies
    Can probably be used for automation
    The official commandline tools may have a better performance
"""


def downloadDataSetbyTaxonomyName(
    taxonomy_name_list: list,
    root_dir_path: str = os.path.join(
        os.getcwd(),
        "database",
    ),
) -> None:
    try:
        for item in taxonomy_name_list:
            assert isinstance(item, str)
    except BaseException:
        raise ValueError(f"Found zero valid taxonomy name in query {', '.join(taxonomy_name_list)}")

    print("Querying taxonomy id")
    genome_api = GenomeApi()
    matched_taxnonmy_name = {}
    request = []
    for taxonomy_name in taxonomy_name_list:
        genome_respond = genome_api.genome_tax_name_query(taxonomy_name, async_req=True)
        request.append(genome_respond)
    for respond in request:
        sci_name_and_ids = respond.get()
        try:
            for entry in sci_name_and_ids["sci_name_and_ids"]:
                matched_taxnonmy_name[entry["tax_id"]] = entry["sci_name"]
        except ApiAttributeError:
            continue
    print("Matched taxonomy id")
    for key, value in matched_taxnonmy_name.items():
        print(key + "\t" + value)
    if len(matched_taxnonmy_name.keys()) == 0:
        raise RuntimeError(f"No valid assembly can be found with the given query: {' ,'.join(taxonomy_name_list)}")
    downloadDataSetbyTaxonomyID(taxonomy_id_list=list(matched_taxnonmy_name.keys()), root_dir_path=root_dir_path)


def downloadDataSetbyTaxonomyID(taxonomy_id_list: list, root_dir_path: str = os.path.join(os.getcwd(), "database")):
    try:
        for item in taxonomy_id_list:
            assert isinstance(item, str) or isinstance(item, int)
        taxonomy_id_list = [int(item) for item in taxonomy_id_list]
        taxonomy_id_list = [str(item) for item in taxonomy_id_list]
    except BaseException:
        raise ValueError("Input taxonomyID Invalid")

    print("Querying assembly descriptors")
    genome_api = GenomeApi()
    request = []
    assembly_accessions = {}

    for taxonomy_id in taxonomy_id_list:
        genome_respond = genome_api.assembly_descriptors_by_taxon(taxonomy_id, async_req=True, page_size=1000)
        request.append(genome_respond)

    for respond in request:
        assembly_descriptors = respond.get()
        try:
            for item in assembly_descriptors["assemblies"]:
                assembly_accessions[item["assembly"]["assembly_accession"]] = [item["assembly"]["biosample"]["description"]["organism"]["tax_id"], item["assembly"]["biosample"]["description"]["organism"]["organism_name"]]
        except ApiAttributeError:
            continue
    print("Matched assembly accession")
    for key, value in assembly_accessions.items():
        print(f"{str(key)}\t{str(value[0])}\t{str(value[1])}")
    if len(assembly_accessions.values()) == 0:
        raise RuntimeError(f"No valid assembly can be found with the given query: {' ,'.join(taxonomy_id_list)}")
    downloadDataSetbyAccession(list(assembly_accessions.keys()), root_dir_path=root_dir_path)


def downloadDataSetbyAccession(accession_list: list, root_dir_path: str = os.path.join(os.getcwd(), "database")):
    try:
        for item in accession_list:
            assert isinstance(item, str)
    except BaseException:
        raise ValueError(f"Found zero valid assembly accession in query {', '.join(accession_list)}")
    try:
        if not os.path.exists(root_dir_path):
            os.makedirs(root_dir_path)
        elif len(os.listdir(root_dir_path)) > 0:
            print("Removing previous file")
            shutil.rmtree(root_dir_path)
            os.makedirs(root_dir_path)

    except BaseException:
        raise ValueError("root dir path is not valid, please check input")
    genome_api = GenomeApi()
    print("Downloading assemblies")
    downloaded = genome_api.download_assembly_package_post(
        v1_assembly_dataset_request=V1AssemblyDatasetRequest(
            assembly_accessions=list(accession_list),
            include_annotation=True,
            include_annotation_type=[V1AnnotationForAssemblyType("DEFAULT"), V1AnnotationForAssemblyType("GENOME_GFF"), V1AnnotationForAssemblyType("PROT_FASTA")],
            include_tsv=True,
        )
    )
    downloaded.close()
    filename = os.path.split(downloaded.name)[-1]
    filepath = os.path.join(root_dir_path, filename)
    print("Moving downloadedd files to storage path")
    if os.path.exists(filepath):
        os.remove(filepath)
    shutil.move(downloaded.name, filepath)
    # download_respond=genome_api.download_assembly_package_post(list(assembly_accessions.keys()))
    subprocess.run(["unzip", "-d", root_dir_path, "-qo", filepath])
    for item in os.listdir(os.path.join(root_dir_path, os.path.splitext(filename)[0], "data")):
        item_path = os.path.join(root_dir_path, os.path.splitext(filename)[0], "data", item)
        shutil.move(item_path, root_dir_path)
    os.remove(os.path.join(root_dir_path, filename))
    shutil.rmtree(os.path.join(root_dir_path, os.path.splitext(filename)[0]))


if __name__ == "__main__":
    downloadDataSetbyAccession(
        root_dir_path="/home/hewc/programming/python_project_linux/pytarget/database",
        accession_list=[
            "GCA_004770555.1",
        ],
    )
