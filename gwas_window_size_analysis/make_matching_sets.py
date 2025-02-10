# Purpose of file: To perform an analysis to compare the effect of different window sizes on GWAS results
# Author: Sophie Sigfstead
import random
import subprocess
import os
import pandas as pd 
import matplotlib.pyplot as plt
import shutil
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def main():
    # Step 1. Define a dictionary that houses random track sets to run the single track analysis on. 
    # constants to set up the dictionary
    num_track_sets = 19
    num_tracks_per_set = 33
    track_set_dictionary = {}

    # Brain tracks - to be excluded from the random track sets
    brain_tracks = {
        0, 1, 9, 76, 78, 80, 81, 172, 179, 216, 240, 261, 278,
        319, 326, 338, 355, 370, 403, 411, 421, 458, 462, 469,
        499, 524, 552, 580, 582, 602, 644, 669
    }

    # Create the keys - these will be the IDs of the track sets
    keys = [chr(i) for i in range(ord('a'), ord('a') + num_track_sets)]

    # Create the list of tracks to choose from (excluding brain tracks)
    valid_numbers = [i for i in range(684) if i not in brain_tracks]

    # Ensure no overlap between track sets by sampling without replacement
    used_tracks = set()

    for key in keys:
        # Filter the valid numbers to exclude any previously used tracks
        available_tracks = [num for num in valid_numbers if num not in used_tracks]
        
        # Sample without replacement
        sampled_tracks = random.sample(available_tracks, num_tracks_per_set)
        track_set_dictionary[key] = sampled_tracks
        
        # Update the used tracks set
        used_tracks.update(sampled_tracks)

    for key, value in track_set_dictionary.items():
        # Step 1. Match up the summary statistics with the data from 1000 genomes for each track set. 
        # This will create 19 track sets. 
        command = [
            'python', 
            '../gwas_1_matching.py', 
            '../../GWAS_Data/GCST90277450.tsv', 
            '../../GWAS_Data/1000genomes_as_csv', 
            f"{value}",
            f'./gwas_1_matching_random_sets/id={key}'
        ]
        subprocess.run(command)  

    
if __name__ == "__main__":
    main()



