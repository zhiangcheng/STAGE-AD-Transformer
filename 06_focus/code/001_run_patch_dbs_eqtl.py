# This code replaces the rsid column of the db with the varID
# This helps the db when imported to FOCUS format since our LD
# data uses varID ~ chr<num>_<pos>_<a1>_<a2> to identify variants.
import os
import sqlite3
from sqlite3 import OperationalError
import glob


def create_empty_template(original_db, template_db):
    """Copy table schemas (without data) from the source database to the template database."""
    # Connect to the source database
    src = sqlite3.connect(original_db)
    src.row_factory = sqlite3.Row  # Make SQL statements easier to extract

    # Get the CREATE statements for all tables, excluding system tables
    cursor = src.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';")
    create_statements = [row[0] for row in cursor.fetchall() if row[0] is not None]

    # Get the CREATE statements for all indexes
    cursor = src.execute("SELECT sql FROM sqlite_master WHERE type='index' AND name NOT LIKE 'sqlite_%';")
    index_statements = [row[0] for row in cursor.fetchall() if row[0] is not None]

    src.close()

    # Create a new database and execute the table and index creation statements
    if os.path.exists(template_db):
        os.remove(template_db)  # Overwrite the existing file
    dst = sqlite3.connect(template_db)
    for stmt in create_statements:
        dst.execute(stmt)
    for stmt in index_statements:
        dst.execute(stmt)
    dst.commit()
    dst.close()

# Example: generate a template from any source database
os.makedirs("../input/template_db", exist_ok=True)
original = "/media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_eqtl_dbs_mashr/mashr_Adipose_Visceral_Omentum.db"
create_empty_template(original, "../input/template_db/empty_template_eqtl.db")


og_dbs = glob.glob("/media/desk15/iy2120/TWAS/Breast-cancer-Example/data/prediction_model/gtex_v8_eqtl_dbs_mashr/*.db")
num_og = len(og_dbs)

# con = sqlite3.connect("combine_sqtl.db")
# cur = con.cursor()
# print("Connected to main db")

os.makedirs("../input/patched/eqtl_dbs", exist_ok=True)

num_done = 0
for index, single_db in enumerate(og_dbs):
    target = f"../input/patched/eqtl_dbs/{os.path.basename(single_db)}"
    os.system(f"cp ../input/template_db/empty_template_eqtl.db {target}")
    db_fname = os.path.basename(single_db)
    db_tiss_name = db_fname.split('mashr_')[1].split('.')[0] # Depends on mashr_<tissue_name>.db format

    con = sqlite3.connect(target)
    cur = con.cursor()
    
    cur.execute("ATTACH '{}' as db{};".format(single_db, index))
    # print("Attached {}".format(single_db))

    combine = "INSERT INTO " + "extra" + " SELECT * FROM db{}.".format(index) + "extra"
    # print(combine)
    cur.execute(combine)

    combine = "INSERT INTO " + "weights" + " SELECT gene, varID, varID, ref_allele, eff_allele, weight FROM db{ind}.".format(ind=index) + "weights;"
    print(db_tiss_name)
    print(combine)
    print()
    cur.executescript(combine)

    cur.execute("detach database db{};".format(index))
    cur.close()

    con.close()
