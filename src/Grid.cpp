#include "../include/Grid.h"

std::vector<Triangle> Grid::get_unique_tringles() {
    std::set<std::array<int, 3>> iunique_triangles;
    std::vector<Triangle> utriangles;
    for(auto &t: itriangles){
        std::array<int, 3> temp({t.iv1, t.iv2, t.iv3});
        std::sort(temp.begin(), temp.end()); //TODO: возможно лучше ручками
        iunique_triangles.insert(temp);
    }
    utriangles.reserve(iunique_triangles.size());
    for(auto &t: iunique_triangles){
        utriangles.emplace_back(points[t[0]], points[t[1]], points[t[2]]);
    }
    return utriangles;
}

void Grid::read() {
    data_file >> num_points;
    points.resize(num_points);
    for (int i = 0; i < num_points; ++i) {
        data_file >> points[i].x >> points[i].y >> points[i].z;
    }
    data_file >> num_triangles;
    itriangles.reserve(2 * num_triangles);
    for (int i = 0; i < num_triangles; ++i) {
        size_t v1, v2, v3, v4;
        data_file >> v1 >> v2 >> v3 >> v4;
        --v1, --v2, --v3, --v4;
        if(v3 == v4){
            itriangles.emplace_back(v1, v2, v3);
        } else {
            itriangles.emplace_back(v1, v2, v3);
            itriangles.emplace_back(v2, v3, v4);
        }
    }
    num_triangles = itriangles.size();
}

void Grid::insert_edge(int e1, int e2, int v) {
    if(edges.find({e1, e2}) == edges.end()){
        edges[{e1, e2}] = {v, -1};
    } else {
        edges[{e1, e2}].second = v;
    }
}

void Grid::get_unique_edges() {
    for(auto &t : itriangles){
        std::array<int, 3> e({t.iv1, t.iv2, t.iv3});
        std::sort(e.begin(), e.end()); //TODO: возможно лучше ручками
        insert_edge(e[0], e[1], e[2]);
        insert_edge(e[0], e[2], e[1]);
        insert_edge(e[1], e[2], e[0]);
    }
}
void Grid::print_edges() {
    for (auto [k, v]: edges) {
        std::cout << k.first << " " << k.second << ": " << v.first << " " << v.second << '\n';
    }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> Grid::get_coordinates() {
    std::vector<double> x, y, z;
    for(auto p : points){
        x.push_back(p.x);
        y.push_back(p.y);
        z.push_back(p.z);
    }
    return std::make_tuple(x, y, z);
}

Grid::Grid(const char* file_name) : data_file(file_name){
    if(!data_file.is_open()){
        std::cout << "can't open: " << file_name << std::endl;
        return;
    }
}

std::array<vec3, 3> Grid::get_points(iTriangle t) {
    return {vec3(points[t.iv1]), vec3(points[t.iv2]), vec3(points[t.iv3])};
}

std::array<double, 3> Grid::get_point_coord(int i) {
    return {points[i].x, points[i].y, points[i].z};
}

void Grid::getMarkedTriangles() {
    for(auto [e, v] : edges){
        triangles.insert(std::make_pair<std::tuple<int, int, int>, Triangle>
                                 ({e.first, e.second, v.first}, {points[e.first], points[e.second], points[v.first]})
        );
        triangles.insert(std::make_pair<std::tuple<int, int, int>, Triangle>
                                 ({e.first, e.second, v.second}, {points[e.first], points[e.second], points[v.second]})
        );
    }
}

bool Grid::check_dist(pair<int, int> p1, pair<int, int> p2) const { //TODO: в треугольнике можно хранить max ребро
    vec3 m1 = (points[p1.first] + points[p1.second]) / 2;
    vec3 m2 = (points[p2.first] + points[p2.second]) / 2;
    auto v1 = edges.find(p1)->second;
    auto v2 =  edges.find(p2)->second;
    double d1 = max_edge({p1.first, p1.second, v1.first});
    double d2 = max_edge({p1.first, p1.second, v1.second});
    double d3 = max_edge({p2.first, p2.second, v2.first});
    double d4 = max_edge({p2.first, p2.second, v2.second});
    double d = max(max(d1, d2), max(d3, d4));
    return dist(m1, m2) > 2 * d;
}

double Grid::diametr_grid() {
    double ans = -1;
    for(auto &t : itriangles){
        double temp = max_edge(t);
        if(temp > ans) {
            ans = temp;
        }
    }
    return ans;
}
