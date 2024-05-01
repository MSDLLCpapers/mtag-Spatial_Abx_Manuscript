# This changes the config files to name every study name, and then runs plot_heatmap_3.py

import re
import os

# load constants from a file to import WHICH_HEATMAP, the only constant used from bug_config

exec(open('./bug_config.py').read(), globals())

listFileNames = [
    "3402_CECUM", "3402_COLON", "3405_CECUM", "3405_COLON", "3406_CECUM", "3406_COLON", "3407_CECUM", "3407_COLON", "3409_CECUM", "3409_COLON", "3412_CECUM", "3412_COLON", "3413_CECUM", "3413_COLON", "3414_CECUM", "3414_COLON", "3415_CECUM", "3415_COLON", "3422_CECUM", "3422_COLON", "3423_CECUM", "3423_COLON", "3424_CECUM", "3424_COLON", "3428_CECUM", "3428_COLON", "3429_CECUM", "3429_COLON", "3430_CECUM", "3430_COLON", "3434_CECUM", "3434_COLON", "3435_CECUM", "3435_COLON", "3436_CECUM", "3436_COLON", "3438_CECUM", "3438_COLON", "3441_CECUM", "3441_COLON", "3442_CECUM", "3442_COLON", "3446_CECUM", "3446_COLON", "3447_CECUM", "3447_COLON", "3448_CECUM", "3448_COLON", "3451_CECUM", "3451_COLON", "3452_CECUM", "3452_COLON", "3454_CECUM", "3454_COLON"
]
configFilePath = "./bug_config.py"

for fileName in listFileNames:
    print("processing " + fileName)
    f = open(configFilePath, "r+")
    contents = f.read()
    x = re.sub(r"PROCESS_ONE_STUDY\s=\s.*",
               "PROCESS_ONE_STUDY = '" + fileName + "'", contents)
    f.seek(0)
    f.write(x)
    f.truncate()
    f.close()
    os.system(WHICH_HEATMAP) # added a variable to the bug_config.py to specify if I want to get the affinity matrix for the whole study or each inidivudal file
