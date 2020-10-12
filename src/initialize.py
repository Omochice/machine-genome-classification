from pathlib import Path

import yaml


def init_settig(setting: Path) -> None:
    project_dir = setting.parent
    print(f"GENERATING SETTING FILE. -> {str(setting)}")
    priority = [
        "NCBI Taxonomy", "NCBI", "GBIF Backbone Taxonomy", "Catalogue of Life",
        "Open Tree of Life Reference Taxonomy",
        "The Interim Register of Marine and Nonmarine Genera"
    ]
    dummy = {
        "email": "",
        "source_csv": "",
        "gbk_destination": str(project_dir / "data" / "gbk"),
        "taxoninfo_destination": str(project_dir / "data" / "taxon"),
        "graph_destination": str(project_dir / "data" / "img"),
        "weight_destination": str(project_dir / "data" / "weight.json"),
        "graph_pix": 192,
        "priority": priority,
        "focus_rank": "class",
        "invalid_creatures": str(project_dir / "data" / "invalid_creatures.json")
    }
    with open(setting, "w") as f:
        yaml.safe_dump(dummy, f)
    print("generating is done. Please rewrite if you need.")


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    setting = project_dir / "setting.yml"
    if not setting.exists():
        init_settig(setting)


if __name__ == "__main__":
    main()