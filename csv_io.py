import pandas as pd

# read sad track csv
sad_tracks_path = "target_dnase_ataq_tracks_labelled.csv"
sad_tracks = pd.read_csv(sad_tracks_path, delimiter=",", header=0,
                         usecols=["target_ids", "target_labels"])

# Remove modality from label
sad_tracks["target_labels"] = [str(lab).split(sep=":", maxsplit=1)[1] for lab in sad_tracks["target_labels"]]
sad_tracks.rename(columns={"target_labels":"description"}, inplace=True)

# write csv
sad_tracks.to_csv("sad_track_descriptions.csv")

track_description = {id:description for id, description in zip(sad_tracks["target_ids"], sad_tracks["description"])}