#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <time.h>
#include <filesystem>
#include <utility>
using namespace std;


class ColoringProblem
{
public:
    int GetRandom(int a, int b)
    {
        static mt19937 generator;
        uniform_int_distribution<int> uniform(a, b);
        return uniform(generator);
    }

    void ReadGraphFile(string filename)
    {   
        ifstream fin(filename);
        if (fin.is_open()) {
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
                    colors.resize(vertices + 1);
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
        else {
            std::cout << "Unable to open " << filename << "!\n";
            exit(1);
        }
    }

    void GreedyGraphColoringReference()
    {
        vector<int> uncolored_vertices(neighbour_sets.size());
        for (size_t i = 0; i < uncolored_vertices.size(); ++i)
            uncolored_vertices[i] = i;

        while (! uncolored_vertices.empty())
        {
            int index = GetRandom(0, uncolored_vertices.size() - 1);
            int vertex = uncolored_vertices[index];
            int color = GetRandom(1, maxcolor);
            for (int neighbour : neighbour_sets[vertex])
            {
                if (color == colors[neighbour])
                {
                    color = ++maxcolor;
                    break;
                }
            }
            colors[vertex] = color;
            // Move the colored vertex to the end and pop it
            swap(uncolored_vertices[uncolored_vertices.size() - 1], uncolored_vertices[index]);
            uncolored_vertices.pop_back();
        }
    } 

    void GreedyGraphColoringBNB()
    {
        // arbitrary modification for Branch-and-bound algorithm with coloring bounds algorithm
        // that turned out better than the previous versions (x_x) 
        vector<int> uncolored_vertices(neighbour_sets.size());
        for (size_t i = 0; i < uncolored_vertices.size(); ++i)
            uncolored_vertices[i] = i;

        for (auto v_itr = uncolored_vertices.rbegin(); v_itr != uncolored_vertices.rend(); ++v_itr)
        {
            int vertex = *v_itr;

            int color = 1;            
            bool found_min_color = false;
            while (!found_min_color) {

                bool need_repeat = false;
                for (int neighbour : neighbour_sets[vertex])
                {
                    if (color == colors[neighbour])
                    {
                        color += 1;
                        maxcolor = std::max(color, maxcolor);
                        need_repeat = true;
                        break;
                    }
                }
                
                if (!need_repeat) {
                    found_min_color = true;
                }
            }

            colors[vertex] = color;
        }
    } 

    void GreedyGraphColoringImproved() 
    {   
        /* 
            Implementing Welsh-Powell algorithm
            http://mrsleblancsmath.pbworks.com/w/file/fetch/46119304/vertex%20coloring%20algorithm.pdf
        */
        using Vertex = int;
        using VertexIndex = int;
        using Color = int;

        // Fill the vertices list        
        vector<Vertex> vertices(neighbour_sets.size());
        std::iota(vertices.begin(), vertices.end(), 0);

        // Sort the vertices in order of descending degree
        std::sort(vertices.begin(), vertices.end(), 
            [&](Vertex lhs, Vertex rhs) {
                return neighbour_sets[lhs].size() > neighbour_sets[rhs].size();
            }
        );

        auto first_colored_index = vertices.size();
        const auto colorVertex = [&](Vertex v, VertexIndex vi, Color c, vector<Vertex>& group) {
            // swap with the last uncolored vertex in the list
            swap(vertices[vi], vertices[first_colored_index - 1]);
            // move the first_colored_index to the vi
            first_colored_index -= 1;
            // color the vertex
            colors[v] = c;
            // add the vertex to the group
            group.push_back(v);
        };

        const auto areConnected = [&](const Vertex lhs, const Vertex rhs) {
            return neighbour_sets[lhs].count(rhs) != 0;
        };

        Color current_color = 1;
        while (first_colored_index != 0) {
            vector<Vertex> current_group;

            // Color the first vertex in the list with the current color,
            // remove the vertex from the list of uncolored vertices
            VertexIndex lhs = vertices.front();
            colorVertex(lhs, 0, current_color, current_group);
            
            if (first_colored_index > 0) {
                // Go down the list and color every vertex that is not connected to 
                // any of the vertex in the 'current_group' with the same color
                // (NB: top-down iteration is important here)
                for (size_t i = first_colored_index - 1; i > 0; i--) {
                    auto rhs = vertices[i];

                    bool is_connected_to_any_in_the_group = 
                        std::any_of(current_group.cbegin(),
                                    current_group.cend(),
                                    [&](Vertex v) {
                                        return areConnected(rhs, v);
                                    }
                    );
                    if (!is_connected_to_any_in_the_group) {

                        colorVertex(rhs, i, current_color, current_group);
                    }
                }
            }

            // Colored all possible vertices with the current color,
            // so should move to the next color now.
            current_color++;
        }

        maxcolor = current_color;
    }

    bool Check()
    {
        for (size_t i = 0; i < neighbour_sets.size(); ++i)
        {
            if (colors[i] == 0)
            {
                cout << "Vertex " << i + 1 << " is not colored\n";
                return false;
            }
            for (int neighbour : neighbour_sets[i])
            {   
                if (neighbour != i) {
                    if (colors[neighbour] == colors[i])
                    {
                        cout << "Neighbour vertices " << i + 1 << ", " << neighbour + 1 <<  " have the same color\n";
                        return false;
                    }
                }
            }
        }
        return true;
    }

    int GetNumberOfColors()
    {
        return maxcolor;
    }
 
    const vector<int>& GetColors()
    {
        return colors;
    }

private:
    vector<int> colors;
    int maxcolor = 1;
    vector<unordered_set<int>> neighbour_sets;
};

int main(int argc, char** argv)
{
    if (argc != 2) {
        cout << "Usage: ./vertex_coloring <INPUT-FILES-DIRECTORY>\n";
        exit(0);
    }

    const std::filesystem::path inputDirectory{argv[1]};
    ofstream fout("vertex_coloring.csv");

    fout << "Instance; Colors; Time (sec)\n";
    cout << "Instance; Colors; Time (sec)\n";

    for (const auto& inputFile : std::filesystem::recursive_directory_iterator{inputDirectory}) 
    {   
        const auto inputFileName = std::filesystem::absolute(inputFile.path()).string();
        ColoringProblem problem;
        problem.ReadGraphFile(inputFileName);
        clock_t start = clock();
        problem.GreedyGraphColoringBNB();
        if (! problem.Check())
        {
            fout << "*** WARNING: incorrect coloring: ***\n";
            cout << "*** WARNING: incorrect coloring: ***\n";
        }
        fout << inputFileName << "; " << problem.GetNumberOfColors() << "; " << double(clock() - start) / CLOCKS_PER_SEC << '\n';
        cout << inputFileName << "; " << problem.GetNumberOfColors() << "; " << double(clock() - start) / CLOCKS_PER_SEC << '\n';
    }

    fout.close();
    return 0;
}