from primerdesigner.database import Database
from primerdesigner.primer_designer import find

Database.build(input_dir="./example_raw_dataset/",output_dir="./example_database/")
db = Database(config_path="./example_database/config.json")
find(
    db=db,
    include=[
        "Cryptococcus gattii"
    ],
    exclude=[
        "Cryptococcus neoformans"
    ],
    pick_probe=True,
    reference_id="GCF_000185945.1",
    workers=2
)




