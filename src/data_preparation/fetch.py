from Bio import SeqIO
import pandas as pd
import os

from . import fetch_gbk, fetch_taxon, classification
import json
import shutil
from seqtools.seq_tools import gbk_utils
import yaml
from pathlib import Path

if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    data_dst = Path(config["destination"])
    gbk_dst = data_dst / "gbk"
    gbk_dst.mkdir(parents=True, exist_ok=True)

    organelle_df = pd.read_csv(config["source_csv"])

    # csvからすでにダウンロード済みのものを除いたlistを得る
    accs = filter(lambda x: not fetch_gbk.is_exist_gbk(x, gbk_dst),
                  map(fetch_gbk.extract_usable_accession, organelle_df["Replicons"]))

    gbk_cliant = fetch_gbk.OrganelleFetchCliant(config["email"])
    print("gbk file fetching...")
    for record in gbk_cliant.fetch(accs):
        accession = record.name
        # print(accession)
        with open(gbk_dst / f"{accession}.gbk", "w") as f:
            SeqIO.write(record, f, "genbank")

    taxon_cliant = fetch_taxon.TaxonFetchCliant(config["email"])
    json_dst = Path(config["destination"]) / "json"
    json_dst.mkdir(parents=True, exist_ok=True)

    print("taxon fetching...")
    for i, taxon_record in enumerate(taxon_cliant.fetch(gbk_dst.glob("*.gbk"))):
        dst = json_dst / f"taxon_tmp_{i}.json"
        with open(dst, "w") as f:
            json.dump(taxon_record, f)
            # Servers may stop connectiong

    print("concating json...")
    taxon_info = {}
    for taxon_tmp in json_dst.glob("taxon_tmp_*.json"):
        with open(taxon_tmp, "r") as f:
            partial = json.load(f)
        os.remove(taxon_tmp)
        taxon_info.update(partial)
    with open(json_dst / "taxon_info.json", "w") as f:
        json.dump(taxon_info, f, indent=2)

    print("constructing csv...")
    csv_dst = data_dst / "csv"
    csv_dst.mkdir(parents=True, exist_ok=True)
    column_names = [
        "accession", "name", "class", "seq_len", "is_mongrel", "have_contig",
        "complete", "shotgun", "chromosome"
    ]
    df = pd.DataFrame(columns=column_names)
    contigs = []

    for gbk in gbk_dst.glob("*.gbk"):
        for record in SeqIO.parse(gbk, "genbank"):
            acc = record.name
            name = record.description
            cl = classification.get_taxon(acc, config["priority"], taxon_info)
            seq_len = len(record.seq)
            if "contig" in record.annotations:
                have_contig = True
                contig = record.annotations["contig"]
                contigs.append(gbk_utils.parse_contig(contig)["accession"])
            else:
                have_contig = False

            row = [acc, name, cl, seq_len, gbk_utils.is_mongrel(name), have_contig]
            for word in column_names[-3:]:
                row.append(word in name)
            df = df.append(pd.Series(row, index=column_names), ignore_index=True)
    df.to_csv(csv_dst / "creatures.csv")

    print("fetch contig gbks...")
    contig_dst = gbk_dst / "contigs"
    contig_dst.mkdir(exist_ok=True)
    for contig_record in gbk_cliant.fetch(contigs):
        acc = contig_record.name
        dst = contig_dst / f"{acc}.gbk"
        with open(dst, "w") as f:
            SeqIO.write(contig_record, f, "genbank")
