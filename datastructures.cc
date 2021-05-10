// Datastructures.cc

#include "datastructures.hh"

#include <random>
#include <QDebug>
#include <cmath>
#include <queue>

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator
typedef std::pair<int, Node*> myNodePair;
typedef std::pair<int, Way*> myWayPair;

std::vector<std::tuple<Coord, WayID, Distance>> getPath(Node* start, Node* end, Node* next, Node* cycleNode = 0) {
    std::vector<std::tuple<Coord, WayID, Distance>> x = {{NO_COORD, NO_WAY, NO_DISTANCE}};

    if (start->coord == end->coord) {
        if (next != 0) {
            x.clear();
            x.push_back({start->coord, next->routeTaken->wayID, end->distance});
        }
        return  x;
    }
    else if (end->prevNode == 0) {
        return  x;
    }
    else {
        x.clear();
        if (next == 0) {
            if (cycleNode != 0) {
                x.push_back({cycleNode->coord, NO_WAY, end->distance});
                auto way = std::find_if(end->neighbours.begin(), end->neighbours.end(), [&cycleNode](std::pair<Node*, Way*> &p) {return p.first == cycleNode;})->second;
                x.push_back({end->coord, way->wayID, end->distance});
            }
            else {
                x.push_back({end->coord, NO_WAY, end->distance});
            }
        }
        else if (next != 0) {
            x.push_back({end->coord, next->routeTaken->wayID, end->distance});
        }
        auto recursion = getPath(start, end->prevNode, end);
        x.insert(x.end(), recursion.begin(), recursion.end());
        return x;
    }
}

void relax(Node* u, Node* v, Way* routeTaken = 0) {
    if (v->prevNode != 0) {
        if (v->distance > u->distance + routeTaken->distance) {
            v->distance = u->distance + routeTaken->distance;
            v->prevNode = u;
            v->routeTaken = routeTaken;
        }
    }
}

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

// Modify the code below to implement the functionality of the class.
// Also remove comments from the parameter names when you implement
// an operation (Commenting out parameter name prevents compiler from
// warning about unused parameters on operations you haven't yet implemented.)

Datastructures::Datastructures()
{
    // Replace this comment with your implementation
}

Datastructures::~Datastructures()
{
    places_.clear(); // O(n)
    placeNameMap_.clear(); // O(n)
    placeCoordMap_.clear(); // O(n)
    placeIDs_.clear(); // O(n)
    areaIDs_.clear(); // O(n)

    areas_.clear(); // O(n)
    superAreas_.clear(); // O(n)
    subAreas_.clear(); // O(n)

    ways_.clear();
    nodes_.clear();
}

int Datastructures::place_count()
{
    int count = places_.size(); // O(1)
    return count;
}

void Datastructures::clear_all()
{
    places_.clear(); // O(n)
    placeNameMap_.clear(); // O(n)
    placeCoordMap_.clear(); // O(n)
    placeIDs_.clear(); // O(n)
    areaIDs_.clear(); // O(n)

    areas_.clear(); // O(n)
    superAreas_.clear(); // O(n)
    subAreas_.clear(); // O(n)

    ways_.clear();
    nodes_.clear();
}

std::vector<PlaceID> Datastructures::all_places()
{
    return placeIDs_;
}

bool Datastructures::add_place(PlaceID id, const Name& name, PlaceType type, Coord xy)
{
    if (places_.find(id) == places_.end()) { // Theta(1), O(n)

        places_.insert({id, Place{id, name, type, xy}}); // O(1)
        placeNameMap_.insert({name, Place{id, name, type, xy}}); // O(log n)
        placeCoordMap_.insert({xy, Place{id, name, type, xy}}); // O(log n)
        placeIDs_.push_back(id); // O(1)

        return true;
    }
    else {
        return false;
    }
}

std::pair<Name, PlaceType> Datastructures::get_place_name_type(PlaceID id)
{
    auto place = places_.find(id); // Theta(1), O(n)
    if (place != places_.end()) {
        return {place->second.name, place->second.type};
    }
    else {
        return {NO_NAME, PlaceType::NO_TYPE};
    }
}

Coord Datastructures::get_place_coord(PlaceID id)
{
    auto it = places_.find(id); // Theta(1), O(n)
    if (it != places_.end()) {
        return it->second.coord;
    }
    else {
        return NO_COORD;
    }
}

bool Datastructures::add_area(AreaID id, const Name &name, std::vector<Coord> coords)
{
    if (areas_.find(id) == areas_.end()) { // Theta(1), O(n)
        areas_.insert({id, Area{id, name, coords}}); // Theta(1), O(n)
        superAreas_.insert({id, -1}); // Theta(1), O(n)
        subAreas_.insert({id, {}}); // Theta(1), O(n)
        areaIDs_.push_back(id); // O(1)
        return true;
    }
    else {
        return false;
    }
}

Name Datastructures::get_area_name(AreaID id)
{
    auto it = areas_.find(id); // Theta(1), O(n)

    if (it != areas_.end()) {
        return it->second.name;
    }
    else {
        return NO_NAME;
    }
}

std::vector<Coord> Datastructures::get_area_coords(AreaID id)
{
    auto it = areas_.find(id);

    if (it != areas_.end()) {
        return it->second.coord;
    }
    else {
        return {NO_COORD};
    }
}

void Datastructures::creation_finished()
{
    // Replace this comment with your implementation
    // NOTE!! It's quite ok to leave this empty, if you don't need operations
    // that are performed after all additions have been done.
}


std::vector<PlaceID> Datastructures::places_alphabetically()
{
    std::vector<PlaceID> placeIDs;

    for (auto const& elem: placeNameMap_) {
        placeIDs.push_back(elem.second.placeID); // O(log n)
    }

    return placeIDs;
}

std::vector<PlaceID> Datastructures::places_coord_order()
{
    std::vector<PlaceID> placeIDs;

    for (auto const& elem: placeCoordMap_) {
        placeIDs.push_back(elem.second.placeID); // O(log n)
    }

    return placeIDs;
}

std::vector<PlaceID> Datastructures::find_places_name(Name const& name)
{
    std::vector<PlaceID> placeIDs;

    for (auto const& elem: places_) {
        if (elem.second.name == name) {
            placeIDs.push_back(elem.first);
        }
    }

    return {placeIDs};
}

std::vector<PlaceID> Datastructures::find_places_type(PlaceType type)
{
    std::vector<PlaceID> placeIDs;

    for (auto const& elem: places_) {
        if (elem.second.type == type) {
            placeIDs.push_back(elem.first);
        }
    }

    return placeIDs;
}

bool Datastructures::change_place_name(PlaceID id, const Name& newname)
{
    auto pos = places_.find(id); // O(1)

    if (pos != places_.end()) {

        // O(n)
        // Only need to change one value so find_if is fine from multimap
        auto placeNamesIter = std::find_if(placeNameMap_.begin(), placeNameMap_.end(), [&id](const std::pair<Name, Place> &p) {return p.second.placeID == id;});
        if (placeNamesIter != placeNameMap_.end()) {
            auto extracted = placeNameMap_.extract(placeNamesIter); // O(log n)
            extracted.key() = newname;
            placeNameMap_.insert(std::move(extracted)); // O(log n)
        }

        pos->second.name = newname;

        // O(n)
        auto placeCoordIter = std::find_if(placeCoordMap_.begin(), placeCoordMap_.end(), [&id](const std::pair<Coord, Place> &p) {return p.second.placeID == id;});
        if (placeCoordIter != placeCoordMap_.end()) {
            placeCoordIter->second.name = newname;
        }

        return true;
    }

    return false;
}

bool Datastructures::change_place_coord(PlaceID id, Coord newcoord)
{
    auto pos = places_.find(id); // Theta(1)

    if (pos != places_.end()) {
        pos->second.coord = newcoord;

        // O(n)
        auto placeNamesIter = std::find_if(placeNameMap_.begin(), placeNameMap_.end(), [&id](const std::pair<Name, Place> &p) {return p.second.placeID == id;});
        if (placeNamesIter != placeNameMap_.end()) {
            placeNamesIter->second.coord = newcoord;
        }

        // O(n)
        auto placeCoordIter = std::find_if(placeCoordMap_.begin(), placeCoordMap_.end(), [&id](const std::pair<Coord, Place> &p) {return p.second.placeID == id;});
        if (placeCoordIter != placeCoordMap_.end()) {
            placeCoordIter->second.coord = newcoord;
            auto handler = placeCoordMap_.extract(placeCoordIter); // O(log n)
            handler.key() = newcoord;
            placeCoordMap_.insert(std::move(handler)); // O(log n)
        }

        return true;
    }

    return false;
}

std::vector<AreaID> Datastructures::all_areas()
{
    return areaIDs_;
}

bool Datastructures::add_subarea_to_area(AreaID id, AreaID parentid)
{
    if (areas_.find(id) != areas_.end() && areas_.find(parentid) != areas_.end()) { // Theta (1), O(n)
        std::unordered_map<AreaID, AreaID>::iterator superIter;
        std::unordered_map<AreaID, std::vector<AreaID>>::iterator subIter;

        superIter = superAreas_.find(id); // Theta(1), O(n)
        if (superIter != superAreas_.end()) {
            if (superIter->second == -1) {
                superIter->second = parentid; // O(1)
            }
            else {
                return false;
            }
        }
        subIter = subAreas_.find(parentid); // Theta(1), O(n)
        subIter->second.push_back(id);
        return true;

    }
    return false;
}

std::vector<AreaID> Datastructures::subarea_in_areas(AreaID id)
{
    std::vector<AreaID> areas;

    std::unordered_map<AreaID, AreaID>::const_iterator superIter;

    AreaID currentID = id;

    while(true) { // O(n)
        superIter = superAreas_.find(currentID); // Theta(1), O(n)

        if (superIter == superAreas_.end() && areas.empty())
        {
            return {NO_AREA};
        }

        else if (superIter != superAreas_.end()) {
            if (superIter->second != -1) {
                areas.push_back(superIter->second);
                currentID = superIter->second;
            } else {
                return areas;
            }
        } else {
            return areas;
        }
    }
    return areas;
}

std::vector<PlaceID> Datastructures::places_closest_to(Coord xy, PlaceType type)
{
    std::vector<PlaceID> places;

    std::multimap<Coord, Place> mapToSearchFrom;

    // O (n)
    if (type != PlaceType::NO_TYPE) {
        std::unordered_map<PlaceID, Place>::const_iterator iter;
        for (iter = places_.begin(); iter != places_.end(); ++iter) {
            if (iter->second.type == type) {
                mapToSearchFrom.insert({iter->second.coord, iter->second});
            }
        }
    }
    else {
        mapToSearchFrom = placeCoordMap_;
    }

    if (mapToSearchFrom.size() > 3) {
        for (int i = 0; i < 3; ++i) {
            // O(n log n)
            auto const it = std::min_element(mapToSearchFrom.begin(), mapToSearchFrom.end(), [&xy](const std::pair<Coord, Place> &p1, const std::pair<Coord, Place> &p2) {
                return pow(p1.first.x-xy.x, 2)+pow(p1.first.y-xy.y, 2) < pow(p2.first.x-xy.x, 2)+pow(p2.first.y-xy.y, 2);
            });
            places.push_back(it->second.placeID);
            mapToSearchFrom.erase(it);
        }}
    else {
        for (auto const& elem: mapToSearchFrom) {
            places.push_back(elem.second.placeID);
        }

    }

    return places;
}

bool Datastructures::remove_place(PlaceID id)
{
    auto it = places_.find(id);
    if (it != places_.end()) {
        auto const placeNameIter = std::find_if(placeNameMap_.begin(), placeNameMap_.end(), [&id](const std::pair<Name, Place> &p) {return p.second.placeID == id;});
        if (placeNameIter != placeNameMap_.end()) {
            placeNameMap_.erase(placeNameIter);
        }

        auto const placeCoordIter = std::find_if(placeCoordMap_.begin(), placeCoordMap_.end(), [&id](const std::pair<Coord, Place> &p) {return p.second.placeID == id;});
        if (placeCoordIter != placeCoordMap_.end()) {
            placeCoordMap_.erase(placeCoordIter);
        }

        auto const vectorIter = std::find(placeIDs_.begin(), placeIDs_.end(), id);
        if (vectorIter != placeIDs_.end()) {
            placeIDs_.erase(vectorIter);
        }

        places_.erase(it);
        placeIDs_.erase(std::remove(placeIDs_.begin(), placeIDs_.end(), id), placeIDs_.end());
        return true;

    } else {
        return false;
    }
}

std::vector<AreaID> Datastructures::all_subareas_in_area(AreaID id)
{
    std::vector<AreaID> areas;

    std::unordered_map<AreaID, std::vector<AreaID>>::iterator subIter;

    subIter = subAreas_.find(id); // Theta(1)

    if (subIter == subAreas_.end() && areas.empty())
    {
        return {NO_AREA};
    }

    // the idea is to recursively find all the subareas while saving them in the process
    // O(n)
    else if (subIter != subAreas_.end()) {
        std::vector<AreaID>::const_iterator vectorIter;
        for (vectorIter = subIter->second.cbegin(); vectorIter != subIter->second.cend(); ++vectorIter) {
            areas.push_back(*vectorIter);
            auto subAreas_of_subArea = all_subareas_in_area(*vectorIter);
            areas.insert(areas.end(), subAreas_of_subArea.begin(), subAreas_of_subArea.end());
        }
        return areas;
    }
    else {
        return areas;
    }
}

AreaID Datastructures::common_area_of_subareas(AreaID id1, AreaID id2)
{
    auto id1ParentAreas = subarea_in_areas(id1); // O(n)
    auto id2ParentAreas = subarea_in_areas(id2); // O(n)

    std::vector<AreaID> smaller;

    AreaID commonID = NO_AREA;

    auto it = id1ParentAreas.rbegin();
    for (; it != id1ParentAreas.rend(); ++it) { // O(n)
        auto findIt = std::find(id2ParentAreas.begin(), id2ParentAreas.end(), *it); // O(n)
        if (findIt != id2ParentAreas.end()) {
            commonID = *it;
        }
    }

    return commonID;
}

std::vector<WayID> Datastructures::all_ways()
{
    std::vector<WayID> wayIDs;
    for (const auto& way: ways_) {
        wayIDs.push_back(way.first);
    }
    return wayIDs;
}

bool Datastructures::add_way(WayID id, std::vector<Coord> coords)
{
    if (ways_.find(id) == ways_.end()) { // Theta(1), O(n)
        int distance = 0;

        std::vector<Coord>::iterator coordIt;
        for (coordIt = coords.begin(); coordIt != coords.end(); ++coordIt) { // O(y) (amount of coords)
            if (std::next(coordIt, 1) != coords.end()) {
                int dis = sqrt(pow((int)(coordIt->x - std::next(coordIt, 1)->x), 2) + pow((int)(coordIt->y - std::next(coordIt, 1)->y), 2));
                distance += dis;
            }
        }

        std::shared_ptr<Way> newWay(new Way({id, coords, coords.front(), coords.back(), distance}));
        ways_.insert({id, newWay}); // Theta(1), O(n)

        auto firstCoordNodeOld = nodes_.find(coords.front()); // Theta(1), O(n)
        auto secondCoordNodeOld = nodes_.find(coords.back()); // Theta(1), O(n)
        std::pair<std::unordered_map<Coord, std::shared_ptr<Node>, CoordHash>::iterator, bool> node1Iterator;
        std::pair<std::unordered_map<Coord, std::shared_ptr<Node>, CoordHash>::iterator, bool> node2Iterator;

        if (firstCoordNodeOld == nodes_.end()) {
            std::shared_ptr<Node> newNode1(new Node({id+"1", {}, 0, 0, coords.front(), INFINITY, white, 0, -1}));
            node1Iterator = nodes_.insert({coords.front(), newNode1}); // Theta(1), O(n)
        }
        if (secondCoordNodeOld == nodes_.end()) {
            std::shared_ptr<Node> newNode2(new Node({id+"2", {}, 0, 0, coords.back(), INFINITY, white, 0, -1}));
            node2Iterator = nodes_.insert({coords.back(), newNode2}); // Theta(1), O(n)
        }

        // Other one is already
        if (firstCoordNodeOld != nodes_.end() && secondCoordNodeOld == nodes_.end()) {
            firstCoordNodeOld->second->neighbours.push_back({node2Iterator.first->second.get(), newWay.get()});
            node2Iterator.first->second->neighbours.push_back({firstCoordNodeOld->second.get(), newWay.get()});
        }
        // Other one is already
        else if (firstCoordNodeOld == nodes_.end() && secondCoordNodeOld != nodes_.end()) {
            node1Iterator.first->second->neighbours.push_back({secondCoordNodeOld->second.get(), newWay.get()});
            secondCoordNodeOld->second->neighbours.push_back({node1Iterator.first->second.get(), newWay.get()});
        }
        // Both are already stored
        else if (firstCoordNodeOld != nodes_.end() && secondCoordNodeOld != nodes_.end()) {
            if (firstCoordNodeOld == secondCoordNodeOld) {
                return true;
            }
            else {
                firstCoordNodeOld->second->neighbours.push_back({secondCoordNodeOld->second.get(), newWay.get()});
                secondCoordNodeOld->second->neighbours.push_back({firstCoordNodeOld->second.get(), newWay.get()});
            }
        }
        // Both are new nodes
        else if (firstCoordNodeOld == nodes_.end() && secondCoordNodeOld == nodes_.end()) {
            node1Iterator.first->second->neighbours.push_back({node2Iterator.first->second.get(), newWay.get()});
            node2Iterator.first->second->neighbours.push_back({node1Iterator.first->second.get(), newWay.get()});
        }

        return true;
    }
    else {
        return false;
    }
}


std::vector<std::pair<WayID, Coord>> Datastructures::ways_from(Coord xy)
{
    std::vector<std::pair<WayID, Coord>> ways;

    auto nodeIt = nodes_.find(xy);

    if (nodeIt != nodes_.end()) {
        for (auto& neighbour: nodeIt->second->neighbours) {
            auto wayToNeighbour = neighbour.second;
            if (wayToNeighbour->startingCoord == xy) {
                ways.push_back({neighbour.second->wayID, wayToNeighbour->endingCoord});
            } else {
                ways.push_back({neighbour.second->wayID, wayToNeighbour->startingCoord});
            }
        }
    }

    return ways;
}

std::vector<Coord> Datastructures::get_way_coords(WayID id)
{
    auto it = ways_.find(id);

    if (it != ways_.end()) {
        return it->second->coords;
    }
    else {
        return {NO_COORD};
    }
}

void Datastructures::clear_ways()
{
    ways_.clear();
    nodes_.clear();
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_any(Coord fromxy, Coord toxy)
{
    auto startIdIt = nodes_.find(fromxy);
    auto endIdIt = nodes_.find(toxy);
    if (startIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }
    if (endIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }

    std::list<Node*> queue;

    // Mark all nodes as not visited
    for (auto& node: nodes_) {
        node.second->colour = white;
        node.second->distance = INFINITY;
        node.second->prevNode = 0;
        node.second->routeTaken = 0;
        node.second->nodeNoInPath = 0;
        node.second->treeId = -1;
    }

    startIdIt->second->colour = gray;
    startIdIt->second->distance = 0;
    queue.push_back(startIdIt->second.get());
    std::vector<std::tuple<Coord, WayID, Distance>> resultsToReverse;

    // BREADTH FIRST SEARCH

    while (queue.size() != 0) { // O(solmu)
        Node* currentNode = queue.front();
        queue.pop_front();

        for (auto& neighbour: currentNode->neighbours) {
            if (neighbour.first->colour == white) {
                neighbour.first->colour = gray;
                int dis = neighbour.second->distance;
                neighbour.first->distance = dis + currentNode->distance;
                neighbour.first->prevNode = currentNode;
                neighbour.first->routeTaken = neighbour.second;

                queue.push_back(neighbour.first);

                if (neighbour.first->coord == toxy) {
                    resultsToReverse = getPath(startIdIt->second.get(), endIdIt->second.get(), nullptr);
                    goto routeFound;
                }
            }
            currentNode->colour = black;
        }
    }
    routeFound:
    std::vector<std::tuple<Coord, WayID, Distance>> results;

    for (auto reverseIt = resultsToReverse.rbegin(); reverseIt != resultsToReverse.rend(); ++reverseIt) {
        results.push_back(*reverseIt);
    }
    return results;
}

bool Datastructures::remove_way(WayID id)
{
    auto wayIt = ways_.find(id); //Theta(1), O(n)

    if (wayIt != ways_.end()) {
        auto startingNode = nodes_.find(wayIt->second->startingCoord); //Theta(1), O(n)
        auto endingNode = nodes_.find(wayIt->second->endingCoord); //Theta(1), O(n)

        // O(n)
        auto startNeighbourIt = std::find_if(startingNode->second->neighbours.begin(), startingNode->second->neighbours.end(), [&id](std::pair<Node*, Way*> x){return x.second->wayID == id; });
        if (startNeighbourIt != startingNode->second->neighbours.end()) {
            startingNode->second->neighbours.erase(startNeighbourIt);
        }

        // O(n)
        auto endNeighbourIt = std::find_if(endingNode->second->neighbours.begin(), endingNode->second->neighbours.end(), [&id](std::pair<Node*, Way*> x){return x.second->wayID == id; });
        if (endNeighbourIt != endingNode->second->neighbours.end()) {
            endingNode->second->neighbours.erase(endNeighbourIt);
        }

        // If node has no other neighbours, we delete the node.
        if (startingNode->second->neighbours.size() == 0) {
            nodes_.erase(startingNode); // Theta(1) since one element is removed
        }
        if (endingNode->second->neighbours.size() == 0) {
            nodes_.erase(endingNode); // Theta(1) since one element is removed
        }

        ways_.erase(wayIt); // Theta(1) since one element is removed
        return true;
    }
    else {
        return false;
    }
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_least_crossroads(Coord fromxy, Coord toxy)
{
    auto startIdIt = nodes_.find(fromxy);
    auto endIdIt = nodes_.find(toxy);
    if (startIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }
    if (endIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }

    std::list<Node*> queue;

    // Mark all nodes as not visited
    for (auto& node: nodes_) {
        node.second->colour = white;
        node.second->distance = INFINITY;
        node.second->prevNode = 0;
        node.second->routeTaken = 0;
        node.second->nodeNoInPath = 0;
        node.second->treeId = -1;
    }

    startIdIt->second->colour = gray;
    startIdIt->second->distance = 0;
    startIdIt->second->nodeNoInPath = 1;
    queue.push_back(startIdIt->second.get());

    // BREADTH FIRST SEARCH

    while (queue.size() != 0) { // O(solmu)
        Node* currentNode = queue.front();
        queue.pop_front();

        for (auto& neighbour: currentNode->neighbours) { // O(kaari)
            if (neighbour.first->colour == white) {
                neighbour.first->colour = gray;

                if (neighbour.first->nodeNoInPath == 0 || neighbour.first->nodeNoInPath > currentNode->nodeNoInPath + 1) {
                    int dis = neighbour.second->distance;
                    neighbour.first->distance = dis + currentNode->distance;
                    neighbour.first->prevNode = currentNode;
                    neighbour.first->nodeNoInPath = currentNode->nodeNoInPath + 1;
                    neighbour.first->routeTaken = neighbour.second;
                }

                queue.push_back(neighbour.first);
            }
        }
        currentNode->colour = black;
    }

    std::vector<std::tuple<Coord, WayID, Distance>> resultsToReverse = getPath(startIdIt->second.get(), endIdIt->second.get(), nullptr);
    std::vector<std::tuple<Coord, WayID, Distance>> results;

    for (auto reverseIt = resultsToReverse.rbegin(); reverseIt != resultsToReverse.rend(); ++reverseIt) {
        results.push_back(*reverseIt);
    }
    return results;
}

std::vector<std::tuple<Coord, WayID> > Datastructures::route_with_cycle(Coord fromxy)
{
    auto startIdIt = nodes_.find(fromxy);
    if (startIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY}};
    }

    std::list<Node*> stack;

    // Mark all nodes as not visited
    for (auto& node: nodes_) {
        node.second->colour = white;
        node.second->distance = INFINITY;
        node.second->prevNode = 0;
        node.second->routeTaken = 0;
        node.second->nodeNoInPath = 0;
        node.second->treeId = -1;
    }

    stack.push_back(startIdIt->second.get());

    std::vector<std::tuple<Coord, WayID, Distance>> resultsToReverse;

    // DEPTH FIRST SEARCH

    while (stack.size() != 0) { // O(solmu)
        Node* currentNode = stack.back();
        stack.pop_back();

        if (currentNode->colour == white) {
            currentNode->colour = gray;
            stack.push_back(currentNode);

            for (auto& neighbour: currentNode->neighbours) { // O(kaari)
                if (neighbour.first->colour == white) {
                    int dis = neighbour.second->distance;
                    neighbour.first->distance = dis + currentNode->distance;
                    neighbour.first->prevNode = currentNode;
                    neighbour.first->routeTaken = neighbour.second;
                    neighbour.first->nodeNoInPath = currentNode->nodeNoInPath + 1;
                    stack.push_back(neighbour.first);
                }
                // Route can't be backwards
                else if (neighbour.first->colour == gray) {
                    if (currentNode->prevNode != 0 && neighbour.first->coord != currentNode->prevNode->coord) {
                        resultsToReverse = getPath(startIdIt->second.get(), currentNode, 0, neighbour.first);
                        goto cycleIsFound;
                    }
                }
            }
        }
        else {
            currentNode->colour = black;
        }
    }
    cycleIsFound:

    std::vector<std::tuple<Coord, WayID>> results = {};

    if (resultsToReverse.size() != 0) {
        results.clear();
        for (auto reverseIt = resultsToReverse.rbegin(); reverseIt != resultsToReverse.rend(); ++reverseIt) {
            results.push_back({std::get<Coord>(*reverseIt), std::get<WayID>(*reverseIt)});
        }
    }
    return results;
}

std::vector<std::tuple<Coord, WayID, Distance> > Datastructures::route_shortest_distance(Coord fromxy, Coord toxy)
{
    auto startIdIt = nodes_.find(fromxy);
    auto endIdIt = nodes_.find(toxy);
    if (startIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }
    if (endIdIt == nodes_.end()) {
        return {{NO_COORD, NO_WAY, NO_DISTANCE}};
    }

    std::priority_queue<myNodePair, std::vector<myNodePair>, std::greater<myNodePair>> priorityQueue;

    // Mark all nodes as not visited
    for (auto& node: nodes_) {
        node.second->colour = white;
        node.second->distance = INFINITY;
        node.second->prevNode = 0;
        node.second->routeTaken = 0;
        node.second->nodeNoInPath = 0;
        node.second->treeId = -1;
    }

    startIdIt->second->colour = gray;
    startIdIt->second->distance = 0;
    startIdIt->second->nodeNoInPath = 1;
    priorityQueue.push({0, startIdIt->second.get()});

    // DIJKSTRA

    while (priorityQueue.size() != 0) { // O(solmu)
        Node* currentNode = priorityQueue.top().second;
        priorityQueue.pop();

        if (currentNode->coord == toxy) {
            break;
        }

        for (auto& neighbour: currentNode->neighbours) { // O(kaari)
            if (neighbour.first->colour == white) {
                neighbour.first->colour = gray;

                if (neighbour.first->prevNode == 0) {
                    int dis = neighbour.second->distance;
                    neighbour.first->distance = dis + currentNode->distance;
                    neighbour.first->prevNode = currentNode;
                    neighbour.first->routeTaken = neighbour.second;
                }

                priorityQueue.push({neighbour.first->distance, neighbour.first});
            }
            relax(currentNode, neighbour.first, neighbour.second);
        }
        currentNode->colour = black;
    }
    std::vector<std::tuple<Coord, WayID, Distance>> resultsToReverse = getPath(startIdIt->second.get(), endIdIt->second.get(), nullptr);
    std::vector<std::tuple<Coord, WayID, Distance>> results;

    for (auto reverseIt = resultsToReverse.rbegin(); reverseIt != resultsToReverse.rend(); ++reverseIt) {
        results.push_back(*reverseIt);
    }
    return results;
}

Distance Datastructures::trim_ways()
{
    int dist = 0;

    std::priority_queue<myWayPair, std::vector<myWayPair>, std::greater<myWayPair>> waysPriorityQueue;
    std::unordered_map<int, std::vector<Node*>> trees;

    // Add all ways to priority queue with the shortest one first.
    for (auto& way: ways_) {
        waysPriorityQueue.push({way.second->distance, way.second.get()}); // O(log n) + O(1) (vector insert)
    }

    // O(n)
    for (auto& node: nodes_) {
        node.second->colour = white;
        node.second->distance = INFINITY;
        node.second->prevNode = 0;
        node.second->routeTaken = 0;
        node.second->nodeNoInPath = 0;
        node.second->treeId = -1; // Parent tree id used to check if the nodes are in the same tree
    }

    // In the algorithm we take the shortest distance way and add it to a tree with random id.
    // The next ones ways we take we check if a nodes of that way is already in the same tree in which case we remove it (trim).
    // If the nodes are in different trees, we combine those trees to one single tree and change the other halfs treeids to the first halfs.

    while(!waysPriorityQueue.empty()) { // O(n)
        Way* way = waysPriorityQueue.top().second; // O(1)
        waysPriorityQueue.pop(); // O(log n) + O(1) (vector)

        auto node1 = nodes_.find(way->startingCoord)->second.get();
        auto node2 = nodes_.find(way->endingCoord)->second.get();

        // if the nodes are the same (end and start coord are equal) or the nodes have already been handled with so that they are in the same tree --> remove it
        if (node1 == node2 || (node1->colour == gray && node2->colour == gray && node1->treeId == node2->treeId && node1->treeId != -1 && node2->treeId != -1)) {
            remove_way(way->wayID);
        }
        else {
            node1->colour = gray;
            node2->colour = gray;

            auto it1 = trees.find(node1->treeId);
            auto it2 = trees.find(node2->treeId);

            // if nodes are not yet in any tree
            if (node1->treeId == -1 && node2->treeId == -1) {
                int treeId = random_in_range(1, 999999999);
                trees.insert({treeId, {node1, node2}});
                node1->treeId = treeId;
                node2->treeId = treeId;
            }
            // if in different trees, combine them
            else if (node1->treeId != -1 && node2->treeId != -1){
                // Insert smaller tree to bigger tree
                if(it1->second.size()>it2->second.size()) {
                    it1->second.reserve(it1->second.size()+it2->second.size());
                    it1->second.insert(it1->second.end(),it2->second.begin(),it2->second.end());
                    for (auto& node: it2->second) {
                        node->treeId = node1->treeId;
                    }
                    trees.erase(it2);
                } else {
                    it2->second.reserve(it2->second.size()+it1->second.size());
                    it2->second.insert(it2->second.end(),it1->second.begin(),it1->second.end());
                    for (auto& node: it1->second) {
                        node->treeId = node2->treeId;
                    }
                    trees.erase(it1);
                }
            }
            // if node1 is a new node but node2 is already in a tree.
            else if (node1->treeId == -1) {
                it2->second.push_back(node1);
                node1->treeId = it2->first;
            }
            // if node2 is a new node but node1 is already in a tree
            else if (node2->treeId == -1) {
                it1->second.push_back(node2);
                node2->treeId = it1->first;
            }

            dist += way->distance;
        }
    }

    return dist;
}
