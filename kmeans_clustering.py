import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import numpy as np
import import_data_handling_kmeans_clustering as idh
import os


kmeans_kwargs = {
    "init": "k-means++",
    "n_init": 50,
    "max_iter": 500,
    "random_state": 42
}


# Function to scale a list of lists using the StandardScaler
def scale(list_pairwise_distances_core):
    scaler = StandardScaler()
    return scaler.fit_transform(list_pairwise_distances_core)


def elbow_kmeans(scaled_list_pairwise_distances_core, save_folder):
    # Create a list to store the sum of squared errors for each k
    sse = []
    # Loop through k values from 1 to 10
    for k in range(1, 11):
        # Create a KMeans instance with k clusters: model
        model = KMeans(n_clusters=k, **kmeans_kwargs)
        # Fit model to samples
        model.fit(scaled_list_pairwise_distances_core)
        # Append the inertia to the list of inertias
        sse.append(model.inertia_)
    # Plot sse against k
    plt.plot(range(1, 11), sse, '-o')
    plt.xlabel('Number of clusters, k')
    plt.ylabel('Sum of squared distance, sse')
    plt.xticks(range(1, 11))
    plt.savefig(os.path.join(save_folder, 'elbow_plot_reference_genome_assessment.png'))
    plt.close()

    # Find the elbow point
    kn = KneeLocator(range(1, 11), sse, curve="convex", direction="decreasing")
    return kn.elbow


def kmeans_run(k, scaled_list_pairwise_distances_core):
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k, **kmeans_kwargs)
    # Fit model to samples
    model.fit(scaled_list_pairwise_distances_core)
    # Determine the cluster labels of new_points: labels
    cluster_centres = model.cluster_centers_

    return cluster_centres


def determine_sse(coordinates_test_list, centre):
    sse = 0
    for index, coordinate in enumerate(coordinates_test_list):
        sse += (coordinate - centre[index]) ** 2
    return sse


# Function to determine genome closest to the cluster centre
def closest_genome(cluster_centre_np, genome_dict):
    genome_dist = None
    genome_name = None
    for genome, cluster_coordinate in genome_dict.items():
        dist = determine_sse(cluster_coordinate, cluster_centre_np)
        if genome_dist is None or dist < genome_dist:
            genome_dist = dist
            genome_name = genome
    return genome_name


# Find genomes closest to cluster centres
def find_closest_genomes(cluster_centres, genome_dict):
    closest_genomes = []
    for centre in cluster_centres:
        closest_genomes.append(closest_genome(centre, genome_dict))
    return closest_genomes


def return_reference_genomes(file_path, save_folder):
    genome_data = idh.generate_kmeans_data(file_path)
    scalar_transform = scale(genome_data[1])
    np_genome_data = np.array(scalar_transform)
    genome_dict = {}
    for index, genome in enumerate(genome_data[0]):
        genome_dict[genome] = np_genome_data[index]
    k_to_use = elbow_kmeans(np_genome_data, save_folder)
    cluster_centres = kmeans_run(k_to_use, np_genome_data)
    closest_genomes = find_closest_genomes(cluster_centres, genome_dict)

    return closest_genomes



