#include "../include/Grid.h"

std::vector<Triangle> Grid::get_unique_tringles() {
    std::set<std::array<int, 3>> iunique_triangles;
    std::vector<Triangle> utriangles;
    for(auto &t: triangles){
        std::array<int, 3> temp({t.iv1, t.iv2, t.iv3});
        std::sort(temp.begin(), temp.end());
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
        points[i].x /= 50, points[i].y /= 50, points[i].z /= 50;
    }
    data_file >> num_triangles;
    triangles.reserve(2 * num_triangles);
    for (int i = 0; i < num_triangles; ++i) {
        size_t v1, v2, v3, v4;
        data_file >> v1 >> v2 >> v3 >> v4;
        --v1, --v2, --v3, --v4;
        if(v3 == v4){
            triangles.emplace_back(v1, v2, v3);
        } else {
            triangles.emplace_back(v1, v2, v3);
            triangles.emplace_back(v2, v3, v4);
        }
    }
    num_triangles = triangles.size();
}

void Grid::insert_edge(int e1, int e2, int v) {
    if(edges.find({e1, e2}) == edges.end()){
        edges[{e1, e2}] = {v, -1};
    } else {
        edges[{e1, e2}].second = v;
    }
}

void Grid::get_unique_edges() {
    for(auto &t : triangles){
        std::array<int, 3> e({t.iv1, t.iv2, t.iv3});
        std::sort(e.begin(), e.end()); //TODO: возможно лучше ручками
        insert_edge(e[0], e[1], e[2]);
        insert_edge(e[0], e[2], e[1]);
        insert_edge(e[1], e[2], e[0]);
    }
    //triangles.clear(); //TODO: так лучше не делать)
}

/*void Grid::plot_grid() {
    auto [x, y, z] = get_coordinates();
    matplot::scatter3(x,y,z);
    matplot::show();
}*/

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

Grid::Grid(const string &file_name) : data_file(file_name){
    if(!data_file.is_open()){
        std::cout << "can't open: " << file_name << std::endl;
        return;
    }
}

std::array<Point, 3> Grid::get_points(iTriangle t) {
    return {Point(points[t.iv1]), Point(points[t.iv2]), Point(points[t.iv3])};
}

std::array<double, 3> Grid::get_point_coord(int i) {
    return {points[i].x, points[i].y, points[i].z};
}
