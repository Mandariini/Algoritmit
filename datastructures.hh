// Datastructures.hh

#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <functional>
#include <deque>
#include <map>
#include <set>
#include <cmath>
#include <list>
#include <unordered_set>
#include <memory>

// Types for IDs
using PlaceID = long long int;
using AreaID = long long int;
using Name = std::string;
using WayID = std::string;

// Return values for cases where required thing was not found
PlaceID const NO_PLACE = -1;
AreaID const NO_AREA = -1;
WayID const NO_WAY = "!!No way!!";

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();

// Return value for cases where name values were not found
Name const NO_NAME = "!!NO_NAME!!";

// Enumeration for different place types
// !!Note since this is a C++11 "scoped enumeration", you'll have to refer to
// individual values as PlaceType::SHELTER etc.
enum class PlaceType { OTHER=0, FIREPIT, SHELTER, PARKING, PEAK, BAY, AREA, NO_TYPE };

// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};

// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (sqrt(c1.x*c1.x+c1.y*c1.y) < sqrt(c2.x*c2.x+c2.y*c2.y)) { return true; }
    else if (sqrt(c1.x*c1.x+c1.y*c1.y) == sqrt(c2.x*c2.x+c2.y*c2.y)
             && c2.y < c1.y) { return true; }
    else { return false; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Type for a distance (in metres)
using Distance = int;

// Return value for cases where Duration is unknown
Distance const NO_DISTANCE = NO_VALUE;



// This is the class you are supposed to implement

struct Place {
    PlaceID placeID;
    Name name;
    PlaceType type;
    Coord coord;
};

struct Area {
    AreaID areaID;
    Name name;
    std::vector<Coord> coord;
};

enum colour {
    white = 0,
    gray = 1,
    black = 2
};

struct Way;
struct Node;

struct Way {
    WayID wayID;
    std::vector<Coord> coords;
    Coord startingCoord;
    Coord endingCoord;
    int distance;
};

struct Node {
    std::string id;
    std::vector<std::pair<Node*, Way*>> neighbours;
    Node* prevNode;
    Way* routeTaken;
    Coord coord;
    double distance;
    int colour;
    int nodeNoInPath;
    int treeId; // Used for trim ways
};

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    // Estimate of performance: O(1)
    // Short rationale for estimate: unordered_map command size is constant
    int place_count();

    // Estimate of performance: O(n)
    // Short rationale for estimate: Every containers clear command complexity is linear
    void clear_all();

    // Estimate of performance: O(1)
    // Short rationale for estimate: returns the vector placeIDs_.
    std::vector<PlaceID> all_places();

    // Estimate of performance: ≈Theta(log n), O(n)
    // Short rationale for estimate: Find method into unordered_map is in average Theta(1) but worst case is O(n).
    // Other elements have either O(1) or O(log n) complexities so thats why ≈Theta(log n).
    bool add_place(PlaceID id, Name const& name, PlaceType type, Coord xy);

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: Finding from unordered_map is in average Theta(1)
    std::pair<Name, PlaceType> get_place_name_type(PlaceID id);

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: Finding from unordered_map is in average Theta(1)
    Coord get_place_coord(PlaceID id);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance: O(n)
    // Short rationale for estimate: for loop is O(n)
    std::vector<PlaceID> places_alphabetically();

    // Estimate of performance: O(n)
    // Short rationale for estimate: for loop is O(n)
    std::vector<PlaceID> places_coord_order();

    // Estimate of performance: O(n)
    // Short rationale for estimate: for loop in itself is O(n) and accessing from unordered_map is O(1)
    std::vector<PlaceID> find_places_name(Name const& name);

    // Estimate of performance: O(n)
    // Short rationale for estimate: for loop in itself is O(n) and accessing from unordered_map is O(1)
    std::vector<PlaceID> find_places_type(PlaceType type);

    // Estimate of performance: O(n log n)
    // Short rationale for estimate: Because std::find_if is O(n) and predicate accesses map value which is O(log n)
    bool change_place_name(PlaceID id, Name const& newname);

    // Estimate of performance: O(n log n)
    // Short rationale for estimate: Because std::find_if is O(n) and predicate accesses map value which is O(log n)
    bool change_place_coord(PlaceID id, Coord newcoord);

    // We recommend you implement the operations below only after implementing the ones above

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: unordered_map find and insert methods are in average Theta(1), but worst case is O(n)
    bool add_area(AreaID id, Name const& name, std::vector<Coord> coords);

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: unordered_map find is in average Theta(1), but worst case is O(n) and accessing from unordered_map is O(1)
    Name get_area_name(AreaID id);

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: Std::find from unordered_map is in average Theta(1)
    std::vector<Coord> get_area_coords(AreaID id);

    // Estimate of performance: O(1)
    // Short rationale for estimate: returns the vector areaIDs_.
    std::vector<AreaID> all_areas();

    // Estimate of performance: ≈Theta(1), O(n)
    // Short rationale for estimate: We use unordered_maps which have find operations in average ≈Theta(1), O(n)
    bool add_subarea_to_area(AreaID id, AreaID parentid);

    // Estimate of performance: O(n)
    // Short rationale for estimate: while loop is O(n), inside loop is unordered map find which is in average Theta(1) and in rare cases O(n)
    std::vector<AreaID> subarea_in_areas(AreaID id);

    // Non-compulsory operations

    // Estimate of performance:
    // Short rationale for estimate:
    void creation_finished();

    // Estimate of performance: O(n)
    // Short rationale for estimate: Although we use a for loop and inside it we call the same method again recursively to find the subareas' subareas.
    std::vector<AreaID> all_subareas_in_area(AreaID id);

    // Estimate of performance: O (n log n)
    // Short rationale for estimate: We go over a map with min_element which complexity is O(n) and at every pair we access the map key which is O(log n) with pair.first etc.
    std::vector<PlaceID> places_closest_to(Coord xy, PlaceType type);

    // Estimate of performance: O(n log n)
    // Short rationale for estimate: Find_if goes through map so O(n) and maps value to id which is O(log n) (accessing)
    bool remove_place(PlaceID id);

    // Estimate of performance: O(log n)
    // Short rationale for estimate: log n amount of loops * std::find log(n)
    AreaID common_area_of_subareas(AreaID id1, AreaID id2);

    // Phase 2 operations

    // Estimate of performance: Theta(n)
    // Short rationale for estimate: Loop is linear operation
    std::vector<WayID> all_ways();

    // Estimate of performance: Theta(y), O(n + y)
    // Short rationale for estimate: Unordered map find and insert is Theta(1) and O(n). Looping through coords is Theta(y).
    bool add_way(WayID id, std::vector<Coord> coords);

    // Estimate of performance: O(n)
    // Short rationale for estimate: If the node at that coord has a direct way to every other node it is O(n) where n is amount of ways.
    std::vector<std::pair<WayID, Coord>> ways_from(Coord xy);

    // Estimate of performance: O(n)
    // Short rationale for estimate: Because returning a vector is a linear operation and also finding from umap is Theta(1), O(n).
    std::vector<Coord> get_way_coords(WayID id);

    // Estimate of performance: O(n)
    // Short rationale for estimate: Clearing is a linear operation.
    void clear_ways();

    // Estimate of performance: O(V+E) Vertex, Edge
    // Short rationale for estimate: While loop loops at most V times and neighbour for loop loops at most E times.
    // getPath goes through the path from start to end so its O(E). At start also O(V) to reset nodes.
    std::vector<std::tuple<Coord, WayID, Distance>> route_any(Coord fromxy, Coord toxy);

    // Non-compulsory operations

    // Estimate of performance: Theta(1), O(n)
    // Short rationale for estimate: In case given id is not found, umap find is Theta(1), O(n). If found the operation is O(n)
    // since we loop through the neighbours of nodes determined by the starting and ending coord. In this process we delete nodes from their neighbours.
    bool remove_way(WayID id);

    // Estimate of performance: O(V+E) Vertex, Edge
    // Short rationale for estimate: The same as route_any but we only update node parameters if the new way is shorter in node amount than current path.
    std::vector<std::tuple<Coord, WayID, Distance>> route_least_crossroads(Coord fromxy, Coord toxy);

    // Estimate of performance: O(V+E) Vertex, Edge
    // Short rationale for estimate: While loop O(V) amount of Vertices. For loop O(E) since it goes through all Edges at most two times.
    // getPath goes through the path from start to end so its O(E). At start also O(V) to reset nodes.
    std::vector<std::tuple<Coord, WayID>> route_with_cycle(Coord fromxy);

    // Estimate of performance: O((V+E)lg V)
    // Short rationale for estimate: I used Dijkstra-algorithm to find the shortest distance. As soon as the node we are looking from comes out of the
    // priority queue, we can stop the loop to make it a bit faster.
    std::vector<std::tuple<Coord, WayID, Distance>> route_shortest_distance(Coord fromxy, Coord toxy);

    // Estimate of performance: O(n*log n)
    // Short rationale for estimate: Priority queue loops are O(n) and doind pop() is O(log n) so the total is O(n*log n).
    // Combining two vectors is also linear as well as for loop for changing treeID.
    Distance trim_ways();

private:
    std::unordered_map<PlaceID, Place> places_;
    std::multimap<Name, Place> placeNameMap_;
    std::multimap<Coord, Place> placeCoordMap_;

    std::vector<PlaceID> placeIDs_;
    std::vector<AreaID> areaIDs_;

    std::unordered_map<AreaID, Area> areas_;
    // avaimena alue, arvona alue minkä sisällä tutkittava alue on
    std::unordered_map<AreaID, AreaID> superAreas_;
    // avaimena alue, arvona alue, joka on tutkittava alueen sisällä
    std::unordered_map<AreaID, std::vector<AreaID>> subAreas_;

    ///----------------------------------------------------------------------
    /// PHASE2
    ///----------------------------------------------------------------------
    std::unordered_map<WayID, std::shared_ptr<Way>> ways_;
    std::unordered_map<Coord, std::shared_ptr<Node>, CoordHash> nodes_;
};

#endif // DATASTRUCTURES_HH
