#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <filesystem>
#include <iterator>
using namespace std;


class MaxCliqueProblem
{
public:
    static int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        string line;
        int vertices = 0, edges = 0;
        while (getline(fin, line))
        {
            if (line[0] == 'c')
            {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p')
            {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
            }
            else
            {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }
    }

    void FindCliqueLDF(int randomization, int iterations)
    {
        /*
            Implementing 'Largest Degree First` greedy euristic
            sorting vertices in every candidates subgraph  
        */

        using Vertex = int;
        using VertexIndex = int;
        static mt19937 generator;
    
        for (int iteration = 0; iteration < iterations; ++iteration)
        {   
            // the initial clique is empty
            vector<int> clique;

            // fill in initial candidates list - contains all vertices of the graph
            vector<int> candidates(neighbour_sets.size());
            std::iota(candidates.begin(), candidates.end(), 0);
    
            auto candidates_end_iterator = candidates.end();
            while (candidates_end_iterator != candidates.begin())
            {
                // sort the candidates in order of descending degree
                std::sort(candidates.begin(), candidates_end_iterator, 
                    [&](Vertex lhs, Vertex rhs) {
                        return neighbour_sets[lhs].size() > neighbour_sets[rhs].size();
                    }
                );

                // with the chosen level of randomization, select a random candidate 
                int last = std::distance(candidates.begin(), candidates_end_iterator) - 1;
                int rnd = GetRandom(0, min(randomization - 1, last));
                int vertex = candidates[rnd];

                // add the candidate vertex to the clique
                clique.push_back(vertex);

                // remove candidates that are not connected to the added vertex,
                // (NB: real size of the candidates vector does not change)
                candidates_end_iterator = std::remove_if(
                    candidates.begin(), 
                    candidates_end_iterator, 
                    [&](const auto candidate) {
                        return (neighbour_sets[vertex].count(candidate) == 0);
                    }
                );
            }

            if (clique.size() > best_clique.size())
            {
                best_clique = clique;
            }
        }
    }

    void FindCliqueReference(int randomization, int iterations)
    {
        static mt19937 generator;
        for (int iteration = 0; iteration < iterations; ++iteration)
        {
            vector<int> clique;
            vector<int> candidates(neighbour_sets.size());
            for (int i = 0; i < neighbour_sets.size(); ++i)
            {
                candidates[i] = i;
            }

            shuffle(candidates.begin(), candidates.end(), generator);
            while (! candidates.empty())
            {
                int last = candidates.size() - 1;
                int rnd = GetRandom(0, min(randomization - 1, last));
                int vertex = candidates[rnd];
                clique.push_back(vertex);
                for (int c = 0; c < candidates.size(); ++c)
                {
                    int candidate = candidates[c];
                    if (neighbour_sets[vertex].count(candidate) == 0)
                    {
                        // Move the candidate to the end and pop it
                        swap(candidates[c], candidates[last]);
                        candidates.pop_back();
                        --c;
                    }
                }
                shuffle(candidates.begin(), candidates.end(), generator);
            }
            if (clique.size() > best_clique.size())
            {
                best_clique = clique;
            }
        }
    }

    const vector<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        if (unique(best_clique.begin(), best_clique.end()) != best_clique.end())
        {
            cout << "Duplicated vertices in the clique\n";
            return false;
        }
        for (int i : best_clique)
        {
            for (int j : best_clique)
            {
                if (i != j && neighbour_sets[i].count(j) == 0)
                {
                    cout << "Returned subgraph is not a clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
};

int main(int argc, char** argv)
{
    if (argc != 4) {
        cout << "Usage: ./max_clique_greedy <INPUT-FILES-DIRECTORY> <NUM-ITERATIROS> <RANDOMIZATION>\n";
        exit(0);
    }

    const std::filesystem::path inputDirectory{argv[1]};
    int iterations = std::atoi(argv[2]);
    int randomization = std::atoi(argv[3]);
    ofstream fout("max_clique_greedy.csv");

    fout << "File; Clique; Time (sec)\n";

    for (const auto& inputFile : std::filesystem::recursive_directory_iterator{inputDirectory}) 
    {   
        const auto inputFileName = std::filesystem::absolute(inputFile.path()).string();
        MaxCliqueProblem problem;
        problem.ReadGraphFile(inputFileName);
        clock_t start = clock();
        problem.FindCliqueLDF(randomization, iterations);
        if (! problem.Check())
        {
            cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }
        fout << inputFileName << "; " << problem.GetClique().size() << "; " << double(clock() - start) / CLOCKS_PER_SEC << '\n';
        cout << inputFileName << ", result - " << problem.GetClique().size() << ", time - " << double(clock() - start) / CLOCKS_PER_SEC << '\n';
    }
    fout.close();
    return 0;
}