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
#include <stdexcept>
#include <memory>
#include <iomanip>

#include <mutex>
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <tbb/global_control.h>

namespace utility {
    int GetRandom(int a, int b) {
        static std::mt19937 generator;
        std::uniform_int_distribution<int> udist(a, b);
        return udist(generator);
    }
}

class GraphInfo {
public:
    void ReadGraphFile(const std::string& filename) {
        std::ifstream fin(filename);
        std::string line;
        int vertices = 0, edges = 0;
        while (std::getline(fin, line)) {
            if (line[0] == 'c') {
                continue;
            }

            std::stringstream line_input(line);
            char command;
            if (line[0] == 'p') {
                std::string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
                non_neighbours.resize(vertices);
            } else {
                int start, finish;
                line_input >> command >> start >> finish;
                // Edges in DIMACS file can be repeated, but it is not a problem for our sets
                neighbour_sets[start - 1].insert(finish - 1);
                neighbour_sets[finish - 1].insert(start - 1);
            }
        }

        // fill non_neibours for structure for each vertex
        for (int i = 0; i < vertices; ++i) {
            for (int j = 0; j < vertices; ++j) {
                if (neighbour_sets[i].count(j) == 0 && i != j)
                    non_neighbours[i].insert(j);
            }
        }
    }

    std::size_t getNumberOfVertices() const noexcept {
        return neighbour_sets.size();
    }

    const std::unordered_set<int>& getNeighbours(int vertex) const {
        return neighbour_sets[vertex];
    }

    const std::unordered_set<int>& getNonNeighbours(int vertex) const {
        return non_neighbours[vertex];
    }

private:
    std::vector<std::unordered_set<int>> neighbour_sets;
    std::vector<std::unordered_set<int>> non_neighbours;
};

class MaxCliqueTabuSearch;

class BnBSolver;

class MaxCliqueSolution {
public:
    size_t size() const noexcept{
        return solution_border;
    }

    MaxCliqueSolution(std::shared_ptr<GraphInfo> gptr) : gptr{gptr} {
        // allocate Solution - Free - NonFree structure,
        solution_free_nonfree.resize(gptr->getNumberOfVertices());
        // alocate oldness structure fot the soft tabu technique
        vertex_oldness.resize(gptr->getNumberOfVertices(), std::numeric_limits<int>::max());
        // allocate auxiliary (vertex -> index in SFN) structure
        vertex_index_in_sfn.resize(gptr->getNumberOfVertices(), -1);
        for (size_t i = 0; i < solution_free_nonfree.size(); ++i) {
            solution_free_nonfree[i] = i;
            vertex_index_in_sfn[i] = i;
        }

        ClearClique();
    }

    MaxCliqueSolution& operator=(MaxCliqueSolution&& other) {
        if (this != &other) {
            this->gptr = std::move(other.gptr);
            this->solution_free_nonfree = std::move(other.solution_free_nonfree);
            this->vertex_index_in_sfn = std::move(other.vertex_index_in_sfn);
            this->vertex_tightness = std::move(other.vertex_tightness);
            this->vertex_oldness = std::move(other.vertex_oldness);

            this->solution_border = other.solution_border;
            this->free_border = other.free_border;
        }

        return *this;
    };

    MaxCliqueSolution& operator=(const MaxCliqueSolution& other) = default;
    MaxCliqueSolution(const MaxCliqueSolution& other) = default;

    void ClearClique() {
        solution_border = 0;
        vertex_tightness = std::vector(gptr->getNumberOfVertices(), 0);
        free_border = gptr->getNumberOfVertices();
    }

    int RecomputeTightnessForAll() {   
        // Tightness(i) = how much vertices should be removed from the CLIQUE to insert i
        // (i.e number of non-neigbors of the vertex in the clique)

        // for every vertex i
        for (int i = 0; i < gptr->getNumberOfVertices(); i++) {
            // reset tightness value
            vertex_tightness[i] = 0;

            // increment tightness value for every its non-neigbor in the clique
            for (int qi = 0; qi < solution_border; qi++) {
                int clique_vertex = solution_free_nonfree[qi];
                if (gptr->getNonNeighbours(i).count(clique_vertex) != 0) {
                    // found non-neighbor in the clique
                    vertex_tightness[i] += 1;
                }
            }
        }
    }

    int GetTightness(int vertex) {
        return vertex_tightness[vertex];
    }

    void SwapVertices(int vertex, int border) {
        int vertex_at_border = solution_free_nonfree[border];
        std::swap(solution_free_nonfree[vertex_index_in_sfn[vertex]], solution_free_nonfree[border]);
        std::swap(vertex_index_in_sfn[vertex], vertex_index_in_sfn[vertex_at_border]);
    }

    void InsertToClique(int vertex) {
        SwapVertices(vertex, solution_border++);

        // Inserting vertex to the clique -
        // * tightness of all non-neighbors of the vertex should be incremented
        // * non-neighbors from Free should be moved to Non-Free
        for (int i : gptr->getNonNeighbours(vertex)) {
            int i_idx = vertex_index_in_sfn[i];
            bool is_free = (i_idx >= solution_border) && (i_idx < free_border);
            if (is_free && (GetTightness(i) == 0)) {
                --free_border;
                SwapVertices(i, free_border);
            }

            vertex_tightness[i] += 1;
        }

        vertex_oldness[vertex] = 0;
    }

    void RemoveFromClique(int vertex)
    {
        // Removing vertex from the clique -
        // * tightness of all non-neighbors of the vertex should be decremented
        // * 1-tight non-neighbors from Non-Free should be moved to F (Free)
        for (int i : gptr->getNonNeighbours(vertex)) {
            // skip if vertex `i` is not in Non-Free
            int i_idx = vertex_index_in_sfn[i];
            bool is_nonfree = (i_idx >= free_border);
            if (is_nonfree && (GetTightness(i) == 1)) {
                SwapVertices(i, free_border);
                free_border++;
            }

            vertex_tightness[i] -= 1;
        }
        
        --solution_border;
        SwapVertices(vertex, solution_border);

        vertex_oldness[vertex] = 1;
    }

    void forceInsert(int vertex) {
        for (int i : gptr->getNonNeighbours(vertex)) {
            for (int j_idx = 0; j_idx < solution_border; j_idx++) {
                auto j = solution_free_nonfree[j_idx];
                if (i == j) {
                    // found vertex' neighbour that is in the solution,
                    // removing it to push the vertex
                    RemoveFromClique(i);
                }
            }
        }

        if (GetTightness(vertex) == 0) {
            InsertToClique(vertex);
        } else {
            throw std::runtime_error("Expected zero tightness for the vertex!");
        }
    }

    void updateOldness() {
        for (int i = solution_border; i < solution_free_nonfree.size(); i++) {
            if (vertex_oldness[i] != std::numeric_limits<int>::max()) {
                vertex_oldness[i] += 1; // std::numeric_limits<int>::max() means that the vertex
                                        // has never been in the solution
            }
        }
    }

private:
    std::shared_ptr<GraphInfo> gptr;

    std::vector<int> solution_free_nonfree;
    std::vector<int> vertex_index_in_sfn;
    std::vector<int> vertex_tightness;
    std::vector<int> vertex_oldness;
    int solution_border = 0;
    int free_border = 0;

    friend MaxCliqueTabuSearch;
    friend BnBSolver;
};

class MaxCliqueTabuSearch
{
public:
    MaxCliqueTabuSearch(std::shared_ptr<GraphInfo> gptr) : clique{gptr}, best_clique{gptr}, gptr{gptr} {
    }

    MaxCliqueSolution& getSolution() {
        return clique;
    }

    ////////////////////////
    // Fast ILS Algorithm //
    ////////////////////////

    void RunSearch(int starts, int randomization) {
        // Use ILS (Iterated Local Search) metaheuristics,
        // start from a greedy solution
        RunInitialHeuristic(randomization);

        // then repeatedly execute local search with perturbations
        for (int iteration = 0; iteration < starts; iteration++) {
            auto new_solution_clique = MaxCliqueSolution(clique);
            Perturb(new_solution_clique);
            TwoImprovements(new_solution_clique);
            CheckNewSolution(new_solution_clique);

            if (clique.solution_border > best_clique.solution_border) {
                best_clique = clique;
            }
        }
    }

    void Perturb(MaxCliqueSolution& clique_) {
        static const auto chooseK = [&](){
            int alpha = utility::GetRandom(1, 2 * clique_.solution_border);
            if (alpha != 1) {
                return 1;
            } else {
                // TODO: replace with some distribution function?
                static const std::vector<int> kDist = [](){
                    std::vector<int> res(31);
                    for (int i = 0; i < 16; i++) {
                        res[i] = 2;
                    }
                    for (int i = 16; i < 24; i++) {
                        res[i] = 3;
                    }
                    for (int i = 24; i < 28; i++) {
                        res[i] = 4;
                    }
                    for (int i = 28; i < 30; i++) {
                        res[i] = 5;
                    }
                    for (int i = 30; i < 31; i++) {
                        res[i] = 6;
                    }
                    return res;
                }();

                auto randIdx = utility::GetRandom(0, kDist.size() - 1);
                return kDist[randIdx];
            }
        };

        if (clique_.free_border != clique_.solution_free_nonfree.size()) {
            int k = chooseK();

            if (k == 1) {
                // pick a random non-solution vertex and force insert it
                const auto randNonFreeIdx = utility::GetRandom(clique_.free_border, 
                                                            clique_.solution_free_nonfree.size() - 1);
                clique_.forceInsert(clique_.solution_free_nonfree[randNonFreeIdx]);
            } else {
                // choose what vertices to insert

                // start with a random one
                const auto randNonFreeIdx = utility::GetRandom(clique_.free_border, 
                                                            clique_.solution_free_nonfree.size() - 1);
                clique_.forceInsert(clique_.solution_free_nonfree[randNonFreeIdx]);
                k -= 1;

                // pick i-th vertex among the non-solution vertices within distance 2 
                // (neighbors of its neighbours) from the first (i - 1) vertices

                int lastPickedIdx = randNonFreeIdx;
                int attempts_left = 20;
                while (k > 0) {
                    // pick a random neighbour
                    const auto& first_neighbours = gptr->getNeighbours(randNonFreeIdx);
                    const auto first_neighbour_index = utility::GetRandom(0, first_neighbours.size());
                    // pick a random neighbour of the previously picked random neighbour
                    const auto& second_neighbours = gptr->getNeighbours(first_neighbour_index);
                    const auto second_neighbour_index = utility::GetRandom(0, second_neighbours.size());

                    // force-insert if valid
                    const auto idx_in_sfn = clique_.vertex_index_in_sfn[second_neighbour_index];
                    const bool is_non_solution_vertex = (idx_in_sfn > clique_.free_border);
                    if ((second_neighbour_index != 0) && is_non_solution_vertex) {
                        clique_.forceInsert(clique_.solution_free_nonfree[randNonFreeIdx]);
                        k -= 1;
                    }

                    attempts_left -= 1;
                    if (attempts_left == 0) {
                        break;
                    }
                }
            }
        }

        // insert free vertices into the solution until no free remain
        while (clique_.solution_border < clique_.free_border) {
            const auto free_vertex = clique_.solution_free_nonfree[clique_.free_border - 1];
            clique_.InsertToClique(free_vertex);
        }
    }

    void TwoImprovements(MaxCliqueSolution& clique_) {
        // fill in candidates from the solution vertices
        size_t candidates_left = clique_.solution_border;
        std::vector<int> candidates(candidates_left, 0);
        for (int i = 0; i < clique_.solution_border; i++) {
            candidates[i] = clique_.solution_free_nonfree[i];
        }

        const auto remove_from_candidates = [&](int candidate_idx) {
            candidates_left -= 1;
            std::swap(candidates[candidate_idx], candidates[candidates_left]);
        };

        const auto insert_to_candidates = [&](int candidate_vertex) {
            if (candidates_left == candidates.size()) {
                candidates.push_back(candidate_vertex);
            } else {
                candidates[candidates_left] = candidate_vertex;
            }

            candidates_left += 1;
        };

        while (candidates_left != 0) {
            // lhs candidate is a random solution vertex
            const auto lhs_candidate_idx = utility::GetRandom(0, candidates_left - 1);
            const auto lhs_candidate = candidates[lhs_candidate_idx];

            // rhs candidates are 1-tight non-neighbours of the lhs candidate
            // that are not in the solution
            std::vector<int> rhs_candidates;
            for (auto nonneighbour : gptr->getNonNeighbours(lhs_candidate)) {
                bool is_in_solution = (clique_.vertex_index_in_sfn[nonneighbour] < clique_.solution_border);
                if (!is_in_solution && (clique_.GetTightness(nonneighbour) == 1)) {
                    rhs_candidates.push_back(nonneighbour);
                }
            }

            if (rhs_candidates.size() < 2) {
                remove_from_candidates(lhs_candidate_idx);
                continue;
            }

            bool found_pair = false;
            for (int i_idx = 0; i_idx < rhs_candidates.size(); i_idx++) {
                if (found_pair) break;

                for (int j_idx = 0; j_idx < rhs_candidates.size(); j_idx++) {
                    int i = rhs_candidates[i_idx];
                    int j = rhs_candidates[j_idx];

                    if ((i != j) && gptr->getNeighbours(i).count(j) != 0) {
                        // found a pair of adjacent vertices in rhs candidates,
                        // make a 2-improvement -- remove rhs candidate, insert the rhs pair
                        clique_.RemoveFromClique(lhs_candidate);
                        clique_.InsertToClique(i);
                        clique_.InsertToClique(j);

                        // update candidates, part 1
                        remove_from_candidates(lhs_candidate_idx);
                        insert_to_candidates(i);
                        insert_to_candidates(j);

                        found_pair = true;
                        break;
                    }
                }
            }

            if (!found_pair) {
                remove_from_candidates(lhs_candidate_idx);
            }
        }
    }

    void CheckNewSolution(MaxCliqueSolution& new_solution) {
        int sizeDiff_NewCurrent = new_solution.solution_border - clique.solution_border;
        if (sizeDiff_NewCurrent > 0) {
            clique = std::move(new_solution);
        } else {
            int sizeDiff_NewBest = new_solution.solution_border - best_clique.size();
            float change_probability = 1.0 / static_cast<float>(1 + sizeDiff_NewCurrent * sizeDiff_NewBest);

            std::random_device rd;
            std::mt19937 gen(rd());
            std::discrete_distribution<> d({change_probability, 1.0 - change_probability});
            if (d(gen) == 0) {
                clique = std::move(new_solution);
            }
        }
    }

    virtual std::unordered_set<int> GetClique() const {
        std::unordered_set<int> res;
        for (int i = 0; i < best_clique.solution_border; i++) {
            res.insert(best_clique.solution_free_nonfree[i]);
        }

        return res;
    }

    virtual std::size_t GetCliqueSize() const {
        return best_clique.size();
    }

    virtual bool CheckCorrectness() {
        return _CheckCorrectness(best_clique);
    }

protected:
    friend MaxCliqueSolution;
    std::shared_ptr<GraphInfo> gptr;
    MaxCliqueSolution clique;
    MaxCliqueSolution best_clique;

    bool _CheckCorrectness(const MaxCliqueSolution& clique_) const {
        for (int i_idx = 0; i_idx < clique_.solution_border; i_idx++) {
            for (int j_idx = 0; j_idx < clique_.solution_border; j_idx++) {
                int i = clique_.solution_free_nonfree[i_idx];
                int j = clique_.solution_free_nonfree[j_idx];

                if ((i != j) && (gptr->getNeighbours(i).count(j) == 0))
                {
                    std::cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

    void RunInitialHeuristic(int randomization) {
        /*
            Implementing 'Largest Degree First` greedy euristic
            sorting vertices in every candidates subgraph  
        */

        using Vertex = int;
        using VertexIndex = int;
        static std::mt19937 generator;
        
        int greedy_iterations = 10;
        for (int giteration = 0; giteration < greedy_iterations; ++giteration) {   
            clique.ClearClique();

            // fill in initial candidates list - contains all vertices of the graph
            std::vector<int> candidates(gptr->getNumberOfVertices());
            std::iota(candidates.begin(), candidates.end(), 0);
            
            auto candidates_end_iterator = candidates.end();
            while (candidates_end_iterator != candidates.begin()) {
                // sort the candidates in order of descending degree
                std::sort(candidates.begin(), candidates_end_iterator, 
                    [&](Vertex lhs, Vertex rhs) {
                        return gptr->getNeighbours(lhs).size() > gptr->getNeighbours(rhs).size();
                    }
                );

                // with the chosen level of randomization, select a random candidate 
                int last = std::distance(candidates.begin(), candidates_end_iterator) - 1;
                int rnd = utility::GetRandom(0, std::min(randomization - 1, last));
                int vertex = candidates[rnd];

                // add the candidate vertex to the cliques
                clique.InsertToClique(vertex);

                // remove candidates that are not connected to the added vertex,
                // (NB: real size of the candidates vector does not change)
                candidates_end_iterator = std::remove_if(
                    candidates.begin(), 
                    candidates_end_iterator, 
                    [&](const auto candidate) {
                        return (gptr->getNeighbours(vertex).count(candidate) == 0);
                    }
                );
            }

            // Update best clique
            if (clique.solution_border > best_clique.solution_border) {  
                best_clique = clique;
            }
        }
    }
};

class BnBSolverReference : public MaxCliqueTabuSearch {
public:
    BnBSolverReference(std::shared_ptr<GraphInfo> gptr)
        : MaxCliqueTabuSearch{gptr}
    { };

    void RunBnB(int iterations, int randomization) {
        std::vector<int> candidates(gptr->getNumberOfVertices());
        for (size_t i = 0; i < candidates.size(); ++i)
        {
            candidates[i] = i;
        }
        static std::mt19937 generator;
        shuffle(candidates.begin(), candidates.end(), generator);
        BnBRecursion(candidates);
    }

    std::unordered_set<int> BnBGetClique() const {
        return bnb_best_clique;
    }

    std::size_t BnBGetCliqueSize() const {
        return bnb_best_clique.size();
    }

    bool Check()
    {
        for (int i : bnb_best_clique)
        {
            for (int j : bnb_best_clique)
            {
                if (i != j && gptr->getNeighbours(i).count(j) == 0)
                {
                    std::cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    void BnBRecursion(const std::vector<int>& candidates)
    {
        if (candidates.empty())
        {
            if (bnb_clique.size() > bnb_best_clique.size())
            {
                bnb_best_clique = bnb_clique;
            }
            return;
        }

        if (bnb_clique.size() + candidates.size() <= bnb_best_clique.size())
            return;

        for (size_t c = 0; c < candidates.size(); ++c)
        {
            std::vector<int> new_candidates;
            new_candidates.reserve(candidates.size());
            for (size_t i = c + 1; i < candidates.size(); ++i)
            {
                if (gptr->getNeighbours(candidates[c]).count(candidates[i]) != 0)
                    new_candidates.push_back(candidates[i]);
            }
            bnb_clique.insert(candidates[c]);
            BnBRecursion(new_candidates);
            bnb_clique.erase(candidates[c]);
        }
    }

    std::unordered_set<int> bnb_best_clique;
    std::unordered_set<int> bnb_clique;
};

class BnBSolver : public MaxCliqueTabuSearch {
public:
    BnBSolver(std::shared_ptr<GraphInfo> gptr)
        : MaxCliqueTabuSearch{gptr}
    { };

    void RunBnB(int iterations, int randomization) {
        // run initial euristic to get a good lower bound
        RunSearch(iterations, randomization);
        bnb_best_clique = GetClique();
        
        // fill in initial candidates list all vertices of the graph
        std::vector<int> candidates(gptr->getNumberOfVertices());
        for (size_t idx = 0; idx < candidates.size(); ++idx)
        {
            candidates[idx] = idx;
        }

        // sort the candidates in order of increasing degree
        std::sort(candidates.begin(), candidates.end(), 
            [&](int clhs, int crhs) {
                return gptr->getNeighbours(clhs).size() < gptr->getNeighbours(crhs).size();
            }
        );

        // prepare dummy colors vector
        dummy_colors.resize(gptr->getNumberOfVertices(), 0);
        std::iota(dummy_colors.rbegin(), dummy_colors.rend(), 1);

        // prepare colors vector
        colors.resize(gptr->getNumberOfVertices(), 0);

        std::unordered_set<int> bnb_clique;
        BnBRecursion(candidates, bnb_clique);
    }

    std::unordered_set<int> BnBGetClique() const {
        return bnb_best_clique;
    }

    std::size_t BnBGetCliqueSize() const {
        return bnb_best_clique.size();
    }

    bool Check()
    {
        for (int i : bnb_best_clique)
        {
            for (int j : bnb_best_clique)
            {
                if (i != j && gptr->getNeighbours(i).count(j) == 0)
                {
                    std::cout << "Returned subgraph is not clique\n";
                    return false;
                }
            }
        }
        return true;
    }

private:
    void BnBRecursion(const std::vector<int>& candidates, std::unordered_set<int>& bnb_clique)
    {
        if (candidates.empty())
        {
            const std::lock_guard<std::mutex> lock(bnb_best_clique_mutex); 
            if (bnb_clique.size() > bnb_best_clique.size())
            {
                bnb_best_clique = bnb_clique;
            }
            return;
        }

        if (bnb_clique.size() + candidates.size() <= bnb_best_clique.size())
            return;

        tbb::parallel_for(tbb::blocked_range<int>(0, candidates.size()),
            [&](tbb::blocked_range<int> r) {
                std::unordered_set<int> __range_bnb_clique = bnb_clique;

                for (size_t c = r.begin(); c < r.end(); ++c) {
                    std::vector<int> new_candidates;
                    new_candidates.reserve(candidates.size());
                    for (size_t i = c + 1; i < candidates.size(); ++i) {
                        if (gptr->getNeighbours(candidates[c]).count(candidates[i]) != 0)
                        new_candidates.push_back(candidates[i]);
                    }

                    __range_bnb_clique.insert(candidates[c]);
                    BnBRecursion(new_candidates, __range_bnb_clique);
                    __range_bnb_clique.erase(candidates[c]);
                }       
            }
        );
    }

    void GreedyGraphColoringBNB(const std::vector<int>& uncolored_candidates)
    {
        for (auto v_itr = uncolored_candidates.rbegin(); v_itr != uncolored_candidates.rend(); ++v_itr)
        {
            int vertex = *v_itr;

            int color = 1;            
            bool found_min_color = false;
            while (!found_min_color) {

                bool need_repeat = false;
                for (int neighbour : gptr->getNeighbours(vertex))
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

    std::unordered_set<int> bnb_best_clique;
    std::mutex bnb_best_clique_mutex;
    std::vector<int> colors;
    std::vector<int> dummy_colors;
    int maxcolor = 1;
};

int main(int argc, char** argv)
{
    // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);

    if (argc != 4) {
        std::cout << "Usage: ./max_clique_bnb <INPUT-FILES-DIRECTORY> <NUM-ITERATIROS> <RANDOMIZATION>\n";
        exit(0);
    }

    const std::filesystem::path inputDirectory{argv[1]};
    int iterations = std::atoi(argv[2]);
    int randomization = std::atoi(argv[3]);

    std::ofstream fout("clique_bnb.csv");
    fout << "File; Clique; Time (sec)\n";

    for (const auto& inputFile : std::filesystem::recursive_directory_iterator{inputDirectory}) 
    {   
        const auto inputFileName = std::filesystem::absolute(inputFile.path()).string();

        std::cout << "Processing " << inputFileName << std::endl;

        const auto graphPtr = std::make_shared<GraphInfo>();
        graphPtr->ReadGraphFile(inputFileName);
        BnBSolverReference problem(graphPtr);

        tbb::tick_count t0 = tbb::tick_count::now();

        problem.RunBnB(iterations, randomization);
        if (!problem.Check())
        {
            std::cout << "*** WARNING: incorrect clique ***\n";
            fout << "*** WARNING: incorrect clique ***\n";
        }

        tbb::tick_count t1 = tbb::tick_count::now();

        fout << inputFileName << "; " << problem.BnBGetCliqueSize() << "; " << (t1-t0).seconds() << '\n';
        std::cout << inputFileName << ", result - " << problem.BnBGetCliqueSize() << ", time - " << (t1-t0).seconds() << '\n';
    }

    fout.close();
    return 0;
}