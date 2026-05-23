import re
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

def build_sequence_with_features(components):
    seq_parts = []
    features = []
    position = 1

    for token, label in components:
        text = token.strip()

        if text.startswith("(") and text.endswith(")"):
            fragment = text[1:-1]
        else:
            match = re.fullmatch(r"([A-Z])(\d+)", text)
            if match:
                amino_acid, repeat_count = match.groups()
                fragment = "".join(f"VPG{amino_acid}G" for _ in range(int(repeat_count)))
            else:
                fragment = text

        start = position
        end = position + len(fragment) - 1
        seq_parts.append(fragment)
        features.append({"location": (start, end), "label": label})
        position = end + 1

    return {"seq": "".join(seq_parts), "features": features}

color_sets = [
    "#f58a5e",
    "#faac61",
    "#ffef86",
    "#f8d3a9",
    "#b1ff67",
    "#75c6a9",
    "#b7e6d7",
    "#85dae9",
    "#84b0dc",
    "#9eafd2",
    "#c7b0e3",
    "#ff9ccd",
    "#d6b295",
    "#d59687",
    "#b4abac"
]

def sequence_with_features_to_seqrecord(sequence_data, record_id="elp", name="elp", description="ELP construct"):
    record = SeqRecord(Seq(sequence_data["seq"]), id=record_id, name=name, description=description)
    record.annotations["molecule_type"] = "protein"

    feature_ids = {}

    colors = set(color_sets)

    for feature in sequence_data["features"]:
        start, end = feature["location"]

        idx = ""
        color = ""
        
        if feature["label"] in feature_ids:
            idx, color = feature_ids[feature["label"]]
        else:
            idx = len(feature_ids)
            color = random.choice(tuple(colors))

            colors.remove(color)

            if len(colors) == 0:
                colors = color_sets

            feature_ids[feature["label"]] = (idx, color)


        record.features.append(
            SeqFeature(
                FeatureLocation(start - 1, end),
                type="misc_feature",
                qualifiers={
                    "note": [feature["label"]],
                    "label": [feature["label"]],
                    "feature_id": [idx],
                    "color": [color]
                },
            )
        )

    return record


def write_gp_file(sequence_data, name):
    record = sequence_with_features_to_seqrecord(sequence_data, name=name)
    SeqIO.write(record, f"./output/{name}.gp", "gb")

def humanize_seq(components):
    readable = ""
    elp_repeat = False

    for d, _ in components:
        match = re.fullmatch(r"([A-Z])(\d+)", d)
        if match:
            amino_acid, repeat_count = match.groups()

            if not elp_repeat:
                readable += "-"

            elp_repeat = True

            readable += amino_acid * int(repeat_count)
        else:
            elp_repeat = False
            readable += f"-{d}"

    if readable[0] == "-":
        readable = readable[1:]

    return readable