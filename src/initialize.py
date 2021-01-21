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
    template = {
        "email": "<FILL IN>",
        "source_csv": "<FILL IN>",
        "destination": str(project_dir / "data"),
        "graph_pix": 192,
        "priority": priority,
        "focus_rank": "class",
        "use_limit": 5
    }
    with open(setting, "w") as f:
        yaml.safe_dump(template, f)
    print("generating is done. Please rewrite if you need.")

    (project_dir / "data").mkdir(parents=True, exist_ok=True)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    setting = project_dir / "setting.yml"
    if not setting.exists():
        init_settig(setting)


if __name__ == "__main__":
    main()