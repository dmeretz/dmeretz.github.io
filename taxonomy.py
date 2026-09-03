from itertools import combinations
import pandas as pd

# taxonomy paths: root → order → family → genus → species
taxonomy = {
    "Lychee": ["Angiosperms","Sapindales","Sapindaceae","Litchi","Litchi chinensis"],
    "Rambutan": ["Angiosperms","Sapindales","Sapindaceae","Nephelium","Nephelium lappaceum"],
    "Pulasan": ["Angiosperms","Sapindales","Sapindaceae","Nephelium","Nephelium mutabile"],
    "Longan": ["Angiosperms","Sapindales","Sapindaceae","Dimocarpus","Dimocarpus longan"],
    "Mamoncillo": ["Angiosperms","Sapindales","Sapindaceae","Melicoccus","Melicoccus bijugatus"],
    "Matoa": ["Angiosperms","Sapindales","Sapindaceae","Pometia","Pometia pinnata"],
    "Ackee": ["Angiosperms","Sapindales","Sapindaceae","Blighia","Blighia sapida"],
    "Guarana": ["Angiosperms","Sapindales","Sapindaceae","Paullinia","Paullinia cupana"],
    "Soapberry": ["Angiosperms","Sapindales","Sapindaceae","Sapindus","Sapindus spp."],
    "Maple": ["Angiosperms","Sapindales","Sapindaceae","Acer","Acer spp."],
    "HorseChestnut": ["Angiosperms","Sapindales","Sapindaceae","Aesculus","Aesculus spp."],

    "Lansium": ["Angiosperms","Sapindales","Meliaceae","Lansium","Lansium domesticum"],

    "Mangosteen": ["Angiosperms","Malpighiales","Clusiaceae","Garcinia","Garcinia mangostana"],
    "Bacuri": ["Angiosperms","Malpighiales","Clusiaceae","Platonia","Platonia insignis"],
    "Gulupa": ["Angiosperms","Malpighiales","Passifloraceae","Passiflora","Passiflora edulis"],

    "Soursop": ["Angiosperms","Magnoliales","Annonaceae","Annona","Annona muricata"],
    "SugarApple": ["Angiosperms","Magnoliales","Annonaceae","Annona","Annona squamosa"],
    "Atemoya": ["Angiosperms","Magnoliales","Annonaceae","Annona","Annona × atemoya"],

    "Jaboticaba": ["Angiosperms","Myrtales","Myrtaceae","Plinia","Plinia cauliflora"],
    "ZapoteNegro": ["Angiosperms","Ericales","Ebenaceae","Diospyros","Diospyros digyna"],
}


def taxonomic_distance(a, b):
    path_a = taxonomy[a]
    path_b = taxonomy[b]

    # find lowest common ancestor
    lca_depth = 0
    for x, y in zip(path_a, path_b):
        if x == y:
            lca_depth += 1
        else:
            break

    # upward steps from both nodes to LCA
    return (len(path_a) - lca_depth) + (len(path_b) - lca_depth)


labels = list(taxonomy.keys())
matrix = pd.DataFrame(index=labels, columns=labels, dtype=int)

for a, b in combinations(labels, 2):
    d = taxonomic_distance(a, b)
    matrix.loc[a, b] = d
    matrix.loc[b, a] = d

for a in labels:
    matrix.loc[a, a] = 0

print(matrix)


matrix.to_html("pairwise_distances.html")
matrix.to_csv("pairwise_distances.csv")


from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

condensed = squareform(matrix.values)
Z = linkage(condensed, method="average")


from io import StringIO

def to_newick(Z, labels):
    tree = {i: labels[i] for i in range(len(labels))}
    n = len(labels)

    for i, (a, b, _, _) in enumerate(Z):
        tree[n + i] = f"({tree[int(a)]},{tree[int(b)]})"

    return tree[max(tree.keys())] + ";"

newick = to_newick(Z, labels)
print(newick)


matrix.to_html("pairwise_table.html", border=1)
matrix.to_csv("pairwise_table.csv")

with open("tree.newick", "w") as f:
    f.write(newick)
