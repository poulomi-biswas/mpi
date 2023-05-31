// rmpi.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <mpi.h>

#define ROOT_PROCESS 0

struct Edges { 
    int source, destination;
    double weight;

    Edges() : source(0), destination(0), weight(0.0) {}
    Edges(int s, int d, double w)
        : source(s), destination(d), weight(w) {}
};

std::vector<Edges> readGraph(const std::string& filename) {
    std::vector<Edges> e;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int source, destination;
    while (file >> source >> destination) {
        e.emplace_back(source, destination, 1.0);
    }
    return e;
    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == '#')
            continue;

        std::istringstream iss(line);
        int source, destination;
        double weight;
        if (!(iss >> source >> destination >> weight))
            continue;

        e.emplace_back(source, destination, weight);
    }
    return e;
}

std::vector<Edges> findMST(const std::vector<Edges>& e, int numVertices, int numProcesses, int rank) {
    std::vector<Edges> mst;
    std::vector<bool> visited(numVertices, false);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;

    // Start with vertex 0 in process 0
    if (rank == ROOT_PROCESS) {
        visited[0] = true;
        for (const auto& edge : e) {
            if (edge.source == 0) {
                pq.push(std::make_pair(edge.weight, edge.destination));
            }
        }
    }

    while (!pq.empty()) {
        double weight = pq.top().first;
        int destination = pq.top().second;
        pq.pop();

        if (!visited[destination]) {
            visited[destination] = true;
            mst.emplace_back(destination, destination, weight);

            if (rank == ROOT_PROCESS) {
                for (const auto& edge : e) {
                    if (edge.source == destination && !visited[edge.destination]) {
                        pq.push(std::make_pair(edge.weight, edge.destination));
                    }
                }
            }
        }
    }

    return mst;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        if (rank == ROOT_PROCESS) {
            std::cerr << "Error: Missing input file argument." << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string filename = argv[1];

    // Read graph from file
    std::vector<Edges> e = readGraph(filename);

    // Broadcast number of vertices to all processes
    int numVertices = e.back().source + 1;
    MPI_Bcast(&numVertices, 1, MPI_INT, ROOT_PROCESS, MPI_COMM_WORLD);

    // Compute local MST for each process
    std::vector<Edges> localMST = findMST(e, numVertices, numProcesses, rank);

    // Gather local MSTs to process 0
    std::vector<Edges> mst;
    if (rank == ROOT_PROCESS) {
        mst.resize(numVertices - 1);
    }

    MPI_Gather(localMST.data(), localMST.size() * sizeof(Edges), MPI_BYTE,
        mst.data(), localMST.size() * sizeof(Edges), MPI_BYTE,
        ROOT_PROCESS, MPI_COMM_WORLD);

    // Print MST from process 0
    if (rank == ROOT_PROCESS) {
        std::cout << "Minimum Spanning Tree:" << std::endl;
        for (const auto& edge : mst) {
            std::cout << edge.source << " - " << edge.destination << " : " << edge.weight << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
