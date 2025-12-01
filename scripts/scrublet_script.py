import scrublet as scr
import pandas as pd
import numpy as np

#######################invitro#############################
#count_matrix = invitro
counts_matrix = pd.read_csv('/mnt/8TB/users/shameed/shameed/Doublet predictions/mat_invitro.csv', index_col=0)

print(counts_matrix.head())
# Convert to numpy array (required by Scrublet)
counts_array = counts_matrix.values

# Run Scrublet
scrub = scr.Scrublet(counts_array)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Preview results
print(doublet_scores[:10])
print(predicted_doublets[:10])

results = pd.DataFrame({
    'barcode': counts_matrix.index,
    'doublet_score': doublet_scores,
    'predicted_doublet': predicted_doublets
})

results.to_csv('/mnt/8TB/users/shameed/shameed/Doublet predictions/Scrublet_results_invitro.csv', index=False)

##################invivo###########################
import pandas as pd
import numpy as np
import scrublet as scr
counts_matrix = pd.read_csv('/mnt/8TB/users/shameed/shameed/Doublet predictions/mat_invivo.csv', index_col=0)

print(counts_matrix.head())
# Convert to numpy array (required by Scrublet)
counts_array = counts_matrix.values

# Run Scrublet
scrub = scr.Scrublet(counts_array)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Preview results
print(doublet_scores[:10])
print(predicted_doublets[:10])

results = pd.DataFrame({
    'barcode': counts_matrix.index,
    'doublet_score': doublet_scores,
    'predicted_doublet': predicted_doublets
})

results.to_csv('/mnt/8TB/users/shameed/shameed/Doublet predictions/Scrublet_results_invivo.csv', index=False)
