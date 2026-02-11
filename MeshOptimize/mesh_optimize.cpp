// mesh_optimize.cpp — STL mesh smoother and decimator
// Reads ASCII STL, welds vertices, applies Laplacian smoothing,
// decimates via shortest-edge collapse, recomputes normals, writes ASCII STL.
//
// Build (STL only):
//   g++ -O2 -std=c++17 -fopenmp -o mesh_optimize mesh_optimize.cpp
//
// Build with HDF5 output support:
//   g++ -O2 -std=c++17 -fopenmp -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5
//      -o mesh_optimize mesh_optimize.cpp
//
// Without OpenMP (Apple Clang):
//   g++ -O2 -std=c++17 -o mesh_optimize mesh_optimize.cpp
//
// Usage:  ./mesh_optimize input.stl output.[stl|h5mesh] [smooth_iters] [target_ratio]
//
//   smooth_iters  — number of Laplacian smoothing passes (default: 5)
//   target_ratio  — fraction of triangles to keep, 0.0–1.0 (default: 0.5)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef H5_VERS_MAJOR
// HDF5 already included
#else
#if __has_include(<hdf5.h>)
#include <hdf5.h>
#define HAVE_HDF5 1
#else
#define HAVE_HDF5 0
#endif
#endif

// Minimum vertex/triangle count to engage OpenMP (avoid fork overhead on small meshes)
static constexpr int OMP_THRESH = 10000;

// ───────────────────────── Timer ─────────────────────────

struct Timer {
    using Clock = std::chrono::high_resolution_clock;
    const char* label;
    Clock::time_point t0;
    Timer(const char* l) : label(l), t0(Clock::now()) {}
    ~Timer() {
        double dt = std::chrono::duration<double>(Clock::now() - t0).count();
        std::fprintf(stderr, "  %-28s %.4f s\n", label, dt);
    }
};

// ───────────────────────── Vec3 ─────────────────────────

struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    float dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    float length2() const { return dot(*this); }
    float length() const { return std::sqrt(length2()); }
    Vec3 normalized() const {
        float l = length();
        return l > 1e-12f ? (*this) * (1.0f / l) : Vec3(0, 0, 0);
    }
};

// ───────────────────────── Triangle ─────────────────────

struct Tri {
    int v[3];
    bool removed = false;
};

// ───────────────────────── Mesh ─────────────────────────

struct Mesh {
    std::vector<Vec3> verts;
    std::vector<Tri>  tris;
};

// ─────────── edge key helper (used throughout) ──────────

static inline uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// ──────────────────── Read ASCII STL ────────────────────
// Reads entire file into memory, scans for "vertex " tokens,
// parses coordinates with strtof. ~12x faster than getline+istringstream.

static void parse_stl_raw(const char* path,
                           std::vector<Vec3>& raw_verts,
                           int& ntri) {
    Timer t("parse");
    FILE* f = std::fopen(path, "r");
    if (!f) {
        std::fprintf(stderr, "Error: cannot open %s\n", path);
        std::exit(1);
    }
    std::fseek(f, 0, SEEK_END);
    long sz = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    char* buf = static_cast<char*>(std::malloc(sz + 1));
    std::fread(buf, 1, sz, f);
    buf[sz] = '\0';
    std::fclose(f);

    // Estimate: ~200 bytes per facet → reserve
    raw_verts.reserve(sz / 60);

    char* p = buf;
    char* end = buf + sz;
    while (p < end) {
        p = static_cast<char*>(std::memchr(p, 'v', end - p));
        if (!p) break;
        // check "vertex " (7 chars)
        if (end - p >= 7 && p[1] == 'e' && p[2] == 'r' && p[3] == 't'
            && p[4] == 'e' && p[5] == 'x' && p[6] == ' ') {
            p += 7;
            char* next;
            float x = std::strtof(p, &next); p = next;
            float y = std::strtof(p, &next); p = next;
            float z = std::strtof(p, &next); p = next;
            raw_verts.push_back({x, y, z});
        } else {
            p++;
        }
    }
    ntri = static_cast<int>(raw_verts.size()) / 3;
    std::free(buf);
}

// ──────────────────── Vertex welding ────────────────────
// Quantize coordinates to an integer grid, use uint64 hash map.
// ~1100x faster than the float-key unordered_map approach.

static Mesh weld_vertices(const std::vector<Vec3>& raw_verts, int ntri) {
    Timer t("weld");
    Mesh m;
    const float inv_eps = 1.0f / 1e-6f;

    // Hash: quantize each coord to integer, mix with Knuth multiplicative hash
    auto quantize = [inv_eps](const Vec3& v) -> uint64_t {
        auto qi = [](float val, float ie) -> uint32_t {
            return static_cast<uint32_t>(static_cast<int32_t>(
                std::round(val * ie)));
        };
        uint64_t h = qi(v.x, inv_eps);
        h = h * 2654435761ULL ^ qi(v.y, inv_eps);
        h = h * 2654435761ULL ^ qi(v.z, inv_eps);
        return h;
    };

    std::unordered_map<uint64_t, int> vmap;
    vmap.reserve(raw_verts.size() / 3); // ~1/3 are unique typically
    m.verts.reserve(raw_verts.size() / 3);

    // Remap table: raw vertex index → welded vertex index
    std::vector<int> remap(raw_verts.size());
    for (size_t i = 0; i < raw_verts.size(); i++) {
        uint64_t key = quantize(raw_verts[i]);
        auto [it, inserted] = vmap.try_emplace(key, static_cast<int>(m.verts.size()));
        if (inserted) {
            m.verts.push_back(raw_verts[i]);
        }
        remap[i] = it->second;
    }

    // Build triangles (every 3 consecutive raw verts = 1 triangle)
    m.tris.reserve(ntri);
    for (int i = 0; i < ntri; i++) {
        int base = i * 3;
        int a = remap[base], b = remap[base + 1], c = remap[base + 2];
        if (a == b || b == c || a == c) continue; // degenerate
        m.tris.push_back({{a, b, c}, false});
    }
    return m;
}

// ──────────────────── Read (combined) ───────────────────

static Mesh read_stl(const char* path) {
    std::vector<Vec3> raw_verts;
    int ntri = 0;
    parse_stl_raw(path, raw_verts, ntri);
    std::fprintf(stderr, "  Read %d triangles, %zu raw vertices\n",
                 ntri, raw_verts.size());

    Mesh m = weld_vertices(raw_verts, ntri);
    std::fprintf(stderr, "  After welding: %zu unique vertices, %zu triangles\n",
                 m.verts.size(), m.tris.size());
    return m;
}

// ──────────────────── Write ASCII STL ───────────────────
// Formats triangles in parallel (normal computation + snprintf),
// writes the buffered output sequentially.

static void write_stl(const char* path, const Mesh& m,
                       const char* name = "optimized") {
    Timer t("write");

    // Collect active triangle indices
    std::vector<int> active;
    active.reserve(m.tris.size());
    for (int i = 0; i < static_cast<int>(m.tris.size()); i++)
        if (!m.tris[i].removed) active.push_back(i);
    int n_active = static_cast<int>(active.size());

    // Parallel: compute normals and format each triangle into a buffer
    struct TriBuf { char data[256]; int len; };
    std::vector<TriBuf> bufs(n_active);

    #pragma omp parallel for schedule(static) if(n_active > OMP_THRESH)
    for (int i = 0; i < n_active; i++) {
        auto& tri = m.tris[active[i]];
        const Vec3& a = m.verts[tri.v[0]];
        const Vec3& b = m.verts[tri.v[1]];
        const Vec3& c = m.verts[tri.v[2]];
        Vec3 n = (b - a).cross(c - a).normalized();
        bufs[i].len = std::snprintf(bufs[i].data, sizeof(bufs[i].data),
            "  facet normal %g %g %g\n"
            "     outer loop\n"
            "       vertex %g %g %g\n"
            "       vertex %g %g %g\n"
            "       vertex %g %g %g\n"
            "     endloop\n"
            "  endfacet\n",
            n.x, n.y, n.z,
            a.x, a.y, a.z,
            b.x, b.y, b.z,
            c.x, c.y, c.z);
    }

    // Sequential: write to file
    FILE* f = std::fopen(path, "w");
    if (!f) {
        std::fprintf(stderr, "Error: cannot write %s\n", path);
        std::exit(1);
    }
    std::fprintf(f, "solid %s\n", name);
    for (int i = 0; i < n_active; i++)
        std::fwrite(bufs[i].data, 1, bufs[i].len, f);
    std::fprintf(f, "endsolid %s\n", name);
    std::fclose(f);
    std::fprintf(stderr, "  Wrote %d triangles to %s\n", n_active, path);
}

// ──────────────────── Write HDF5 Mesh ───────────────────
// Writes .h5mesh format: /vertices [N,3], /triangles [M,3], /normals [M,3]
// plus root attributes for format detection.

#if HAVE_HDF5
static void write_h5mesh(const char* path, const Mesh& m) {
    Timer t("write h5mesh");

    // Collect active triangles and reindex vertices
    std::vector<int> active;
    active.reserve(m.tris.size());
    for (int i = 0; i < static_cast<int>(m.tris.size()); i++)
        if (!m.tris[i].removed) active.push_back(i);
    int n_active = static_cast<int>(active.size());

    // Find which vertices are actually used, build old→new index map
    std::vector<int> vert_map(m.verts.size(), -1);
    int n_verts = 0;
    for (int i = 0; i < n_active; i++) {
        auto& tri = m.tris[active[i]];
        for (int k = 0; k < 3; k++) {
            if (vert_map[tri.v[k]] < 0)
                vert_map[tri.v[k]] = n_verts++;
        }
    }

    // Build flat vertex array [N, 3]
    std::vector<float> verts_flat(n_verts * 3);
    for (int i = 0; i < static_cast<int>(m.verts.size()); i++) {
        if (vert_map[i] >= 0) {
            int idx = vert_map[i] * 3;
            verts_flat[idx]     = m.verts[i].x;
            verts_flat[idx + 1] = m.verts[i].y;
            verts_flat[idx + 2] = m.verts[i].z;
        }
    }

    // Build flat triangle index array [M, 3] and normals [M, 3]
    std::vector<int32_t> tris_flat(n_active * 3);
    std::vector<float> normals_flat(n_active * 3);

    #pragma omp parallel for schedule(static) if(n_active > OMP_THRESH)
    for (int i = 0; i < n_active; i++) {
        auto& tri = m.tris[active[i]];
        int a = vert_map[tri.v[0]];
        int b = vert_map[tri.v[1]];
        int c = vert_map[tri.v[2]];
        tris_flat[i * 3]     = a;
        tris_flat[i * 3 + 1] = b;
        tris_flat[i * 3 + 2] = c;

        const Vec3& va = m.verts[tri.v[0]];
        const Vec3& vb = m.verts[tri.v[1]];
        const Vec3& vc = m.verts[tri.v[2]];
        Vec3 n = (vb - va).cross(vc - va).normalized();
        normals_flat[i * 3]     = n.x;
        normals_flat[i * 3 + 1] = n.y;
        normals_flat[i * 3 + 2] = n.z;
    }

    // Create HDF5 file
    hid_t file_id = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::fprintf(stderr, "Error: cannot create HDF5 file %s\n", path);
        std::exit(1);
    }

    // Write /vertices dataset [N, 3] float32
    {
        hsize_t dims[2] = {static_cast<hsize_t>(n_verts), 3};
        hid_t space = H5Screate_simple(2, dims, NULL);
        hid_t dset = H5Dcreate2(file_id, "vertices", H5T_IEEE_F32LE,
                                 space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 verts_flat.data());
        H5Dclose(dset);
        H5Sclose(space);
    }

    // Write /triangles dataset [M, 3] int32
    {
        hsize_t dims[2] = {static_cast<hsize_t>(n_active), 3};
        hid_t space = H5Screate_simple(2, dims, NULL);
        hid_t dset = H5Dcreate2(file_id, "triangles", H5T_STD_I32LE,
                                 space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 tris_flat.data());
        H5Dclose(dset);
        H5Sclose(space);
    }

    // Write /normals dataset [M, 3] float32
    {
        hsize_t dims[2] = {static_cast<hsize_t>(n_active), 3};
        hid_t space = H5Screate_simple(2, dims, NULL);
        hid_t dset = H5Dcreate2(file_id, "normals", H5T_IEEE_F32LE,
                                 space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 normals_flat.data());
        H5Dclose(dset);
        H5Sclose(space);
    }

    // Write root attributes
    {
        // format = "h5mesh"
        hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, 7); // "h5mesh" + null
        H5Tset_strpad(str_type, H5T_STR_NULLTERM);
        hid_t scalar = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(file_id, "format", str_type, scalar,
                                 H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_type, "h5mesh");
        H5Aclose(attr);
        H5Sclose(scalar);
        H5Tclose(str_type);
    }
    {
        // version = 1
        hid_t scalar = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(file_id, "version", H5T_STD_I32LE, scalar,
                                 H5P_DEFAULT, H5P_DEFAULT);
        int ver = 1;
        H5Awrite(attr, H5T_NATIVE_INT, &ver);
        H5Aclose(attr);
        H5Sclose(scalar);
    }
    {
        // num_vertices = N
        hid_t scalar = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(file_id, "num_vertices", H5T_STD_I32LE, scalar,
                                 H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT, &n_verts);
        H5Aclose(attr);
        H5Sclose(scalar);
    }
    {
        // num_triangles = M
        hid_t scalar = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(file_id, "num_triangles", H5T_STD_I32LE, scalar,
                                 H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT, &n_active);
        H5Aclose(attr);
        H5Sclose(scalar);
    }

    H5Fclose(file_id);
    std::fprintf(stderr, "  Wrote %d vertices, %d triangles to %s\n",
                 n_verts, n_active, path);
}
#endif // HAVE_HDF5

// ──────────────── Output format detection ────────────────

static bool has_extension(const char* path, const char* ext) {
    const char* dot = std::strrchr(path, '.');
    if (!dot) return false;
    return std::strcmp(dot + 1, ext) == 0;
}

// ─────────────────── Build adjacency ────────────────────
// Uses vector<vector<int>> with sort+unique per neighbor list.
// Deduplication is parallelized (each vertex list is independent).

using AdjList = std::vector<std::vector<int>>;

static AdjList build_adjacency(const Mesh& m) {
    Timer t("build adjacency");
    int nv = static_cast<int>(m.verts.size());
    AdjList adj(nv);

    // Sequential: collect neighbors (writes to shared adj[v] vectors)
    for (auto& tri : m.tris) {
        if (tri.removed) continue;
        for (int i = 0; i < 3; i++) {
            int a = tri.v[i], b = tri.v[(i + 1) % 3];
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
    }

    // Parallel: sort+unique each vertex's neighbor list independently
    #pragma omp parallel for schedule(dynamic, 256) if(nv > OMP_THRESH)
    for (int i = 0; i < nv; i++) {
        std::sort(adj[i].begin(), adj[i].end());
        adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
    }
    return adj;
}

// ─────────────── Laplacian smoothing ────────────────────
// The per-vertex update reads from m.verts[neighbors] (shared read)
// and writes to new_pos[i] (exclusive write) — embarrassingly parallel.

static void laplacian_smooth(Mesh& m, int iterations, float lambda = 0.5f) {
    AdjList adj = build_adjacency(m);
    int nv = static_cast<int>(m.verts.size());

    // Find boundary vertices (edges with only one adjacent triangle)
    // Using vector<char> instead of vector<bool> for thread-safe reads
    std::vector<char> boundary(nv, 0);
    {
        std::unordered_map<uint64_t, int> ecnt;
        ecnt.reserve(m.tris.size() * 3);
        for (auto& tri : m.tris) {
            if (tri.removed) continue;
            for (int i = 0; i < 3; i++)
                ecnt[edge_key(tri.v[i], tri.v[(i + 1) % 3])]++;
        }
        for (auto& [key, cnt] : ecnt) {
            if (cnt == 1) {
                boundary[static_cast<int>(key >> 32)] = 1;
                boundary[static_cast<int>(key & 0xFFFFFFFF)] = 1;
            }
        }
    }

    std::vector<Vec3> new_pos(nv);
    for (int iter = 0; iter < iterations; iter++) {
        #pragma omp parallel for schedule(static) if(nv > OMP_THRESH)
        for (int i = 0; i < nv; i++) {
            if (boundary[i] || adj[i].empty()) {
                new_pos[i] = m.verts[i];
                continue;
            }
            Vec3 avg{0, 0, 0};
            for (int nb : adj[i])
                avg = avg + m.verts[nb];
            float inv_n = 1.0f / static_cast<float>(adj[i].size());
            Vec3 delta = avg * inv_n - m.verts[i];
            new_pos[i] = m.verts[i] + delta * lambda;
        }
        std::swap(m.verts, new_pos);
    }
}

// ────────────── Edge-collapse decimation ────────────────
// Greedy shortest-edge collapse with vertex versioning.
// Sorted-array scan replaces PQ for cache-friendly access on large meshes.
// New edges from collapses go into a secondary buffer, sorted between passes.

struct EdgeEntry {
    float len2;  // squared length — avoids sqrt, same ordering
    int v0, v1;
    int ver0, ver1;
    bool operator<(const EdgeEntry& o) const { return len2 < o.len2; }
};

static void decimate(Mesh& m, float target_ratio) {
    int active_tris = 0;
    for (auto& t : m.tris)
        if (!t.removed) active_tris++;
    int target = std::max(4, static_cast<int>(active_tris * target_ratio));
    if (active_tris <= target) return;

    Timer timer("decimate");
    std::fprintf(stderr, "  Decimating from %d to ~%d triangles ...\n",
                 active_tris, target);

    int nv = static_cast<int>(m.verts.size());

    // Vertex-to-triangle map
    std::vector<std::vector<int>> v2t(nv);
    for (int i = 0; i < static_cast<int>(m.tris.size()); i++) {
        if (m.tris[i].removed) continue;
        for (int k = 0; k < 3; k++)
            v2t[m.tris[i].v[k]].push_back(i);
    }

    // Per-vertex state
    std::vector<char> v_removed(nv, 0);
    std::vector<int>  v_version(nv, 0);

    // Build initial edge list — deduplicated, lengths computed in parallel
    std::vector<uint64_t> init_keys;
    init_keys.reserve(m.tris.size() * 3);
    for (auto& t : m.tris) {
        if (t.removed) continue;
        init_keys.push_back(edge_key(t.v[0], t.v[1]));
        init_keys.push_back(edge_key(t.v[1], t.v[2]));
        init_keys.push_back(edge_key(t.v[2], t.v[0]));
    }
    std::sort(init_keys.begin(), init_keys.end());
    init_keys.erase(std::unique(init_keys.begin(), init_keys.end()),
                    init_keys.end());

    int ne = static_cast<int>(init_keys.size());
    std::vector<EdgeEntry> edges(ne);

    #pragma omp parallel for schedule(static) if(ne > OMP_THRESH)
    for (int i = 0; i < ne; i++) {
        int a = static_cast<int>(init_keys[i] >> 32);
        int b = static_cast<int>(init_keys[i] & 0xFFFFFFFF);
        edges[i] = {(m.verts[a] - m.verts[b]).length2(), a, b, 0, 0};
    }
    init_keys.clear();
    init_keys.shrink_to_fit();

    // Sort initial edges by length (cache-friendly sequential scan later)
    std::sort(edges.begin(), edges.end());

    // Scratch buffers
    std::vector<int> nbrs0, nbrs1, unique_nbrs;
    std::vector<EdgeEntry> new_edges;
    new_edges.reserve(edges.size() / 2);

    // Multi-pass: scan sorted edges, collapse valid ones, collect new edges.
    // When current pass exhausted, sort new edges and repeat.
    int cursor = 0;
    while (active_tris > target) {
        // Scan current sorted edge list
        bool did_work = false;
        while (cursor < static_cast<int>(edges.size()) && active_tris > target) {
            auto& e = edges[cursor++];
            if (v_removed[e.v0] || v_removed[e.v1]) continue;
            if (v_version[e.v0] != e.ver0 || v_version[e.v1] != e.ver1) continue;

            // Link condition
            nbrs0.clear();
            nbrs1.clear();
            for (int ti : v2t[e.v0]) {
                if (m.tris[ti].removed) continue;
                for (int k = 0; k < 3; k++) {
                    int vv = m.tris[ti].v[k];
                    if (vv != e.v0) nbrs0.push_back(vv);
                }
            }
            for (int ti : v2t[e.v1]) {
                if (m.tris[ti].removed) continue;
                for (int k = 0; k < 3; k++) {
                    int vv = m.tris[ti].v[k];
                    if (vv != e.v1) nbrs1.push_back(vv);
                }
            }
            std::sort(nbrs0.begin(), nbrs0.end());
            nbrs0.erase(std::unique(nbrs0.begin(), nbrs0.end()), nbrs0.end());
            std::sort(nbrs1.begin(), nbrs1.end());
            nbrs1.erase(std::unique(nbrs1.begin(), nbrs1.end()), nbrs1.end());
            int shared = 0;
            {
                auto it0 = nbrs0.begin(), it1 = nbrs1.begin();
                while (it0 != nbrs0.end() && it1 != nbrs1.end()) {
                    if (*it0 == *it1) { shared++; ++it0; ++it1; }
                    else if (*it0 < *it1) ++it0;
                    else ++it1;
                }
            }
            if (shared > 2) continue;

            // Collapse v1 → v0
            int v0 = e.v0, v1 = e.v1;
            m.verts[v0] = (m.verts[v0] + m.verts[v1]) * 0.5f;
            v_version[v0]++;
            v_removed[v1] = 1;

            for (int ti : v2t[v1]) {
                if (m.tris[ti].removed) continue;
                auto& t = m.tris[ti];
                bool has_v0 = (t.v[0] == v0 || t.v[1] == v0 || t.v[2] == v0);
                for (int k = 0; k < 3; k++)
                    if (t.v[k] == v1) t.v[k] = v0;
                if (has_v0) {
                    t.removed = true;
                    active_tris--;
                } else {
                    v2t[v0].push_back(ti);
                }
            }
            v2t[v1].clear();

            // Compact v2t[v0]
            {
                auto& list = v2t[v0];
                int w = 0;
                for (int r = 0; r < static_cast<int>(list.size()); r++)
                    if (!m.tris[list[r]].removed) list[w++] = list[r];
                list.resize(w);
            }

            // Collect new edges
            unique_nbrs.clear();
            for (int ti : v2t[v0]) {
                auto& t = m.tris[ti];
                for (int k = 0; k < 3; k++)
                    if (t.v[k] != v0) unique_nbrs.push_back(t.v[k]);
            }
            std::sort(unique_nbrs.begin(), unique_nbrs.end());
            unique_nbrs.erase(std::unique(unique_nbrs.begin(), unique_nbrs.end()),
                              unique_nbrs.end());
            for (int nb : unique_nbrs) {
                float l2 = (m.verts[v0] - m.verts[nb]).length2();
                new_edges.push_back({l2, v0, nb, v_version[v0], v_version[nb]});
            }
            did_work = true;
        }

        // Current pass exhausted. Switch to new edges for next pass.
        if (new_edges.empty() || !did_work) break;
        edges = std::move(new_edges);
        new_edges.clear();
        new_edges.reserve(edges.size());
        std::sort(edges.begin(), edges.end());
        cursor = 0;
    }

    std::fprintf(stderr, "  Decimation complete: %d triangles remain\n", active_tris);
}

// ─────────────────────── main ───────────────────────────

static void print_help(const char* prog) {
    std::fprintf(stderr,
        "mesh_optimize — smooth and decimate ASCII STL meshes\n"
        "\n"
        "Usage:\n"
        "  %s input.stl output.[stl|h5mesh] [smooth_iters] [target_ratio]\n"
        "  %s -h | --help\n"
        "\n"
        "Positional arguments:\n"
        "  input.stl      Input mesh in ASCII STL format\n"
        "  output         Output mesh path (ASCII STL or HDF5 mesh)\n"
        "                 Extension .h5mesh or .h5 writes HDF5 format;\n"
        "                 any other extension writes ASCII STL.\n"
        "\n"
        "Optional arguments:\n"
        "  smooth_iters   Number of Laplacian smoothing passes (default: 5)\n"
        "                 Each pass moves interior vertices toward the average\n"
        "                 of their neighbors. Boundary vertices are pinned.\n"
        "                 Set to 0 to skip smoothing entirely.\n"
        "\n"
        "  target_ratio   Fraction of triangles to keep after decimation,\n"
        "                 between 0.0 and 1.0 (default: 0.5)\n"
        "                 Uses greedy shortest-edge collapse.\n"
        "                 Set to 1.0 to skip decimation entirely.\n"
        "\n"
        "Processing pipeline:\n"
        "  1. Read    — parse ASCII STL and weld duplicate vertices\n"
        "  2. Smooth  — iterative Laplacian smoothing (preserves boundaries)\n"
        "  3. Decimate — shortest-edge collapse with topological checks\n"
        "  4. Write   — output ASCII STL with recomputed face normals\n"
        "\n"
        "Examples:\n"
        "  %s mesh.stl out.stl              # 5 smooth passes, keep 50%%\n"
        "  %s mesh.stl out.stl 10 0.3       # 10 passes, keep 30%%\n"
        "  %s mesh.stl out.stl 20 1.0       # smooth only, no decimation\n"
        "  %s mesh.stl out.stl 0 0.5        # decimate only, no smoothing\n"
        "  %s mesh.stl out.stl 0 1.0        # weld + recompute normals only\n",
        prog, prog, prog, prog, prog, prog, prog);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help(argv[0]);
        return 1;
    }

    // Check for help flag
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-h") == 0 ||
            std::strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        }
    }

    if (argc < 3) {
        std::fprintf(stderr, "Error: missing output file\n\n");
        print_help(argv[0]);
        return 1;
    }

    const char* infile = argv[1];
    const char* outfile = argv[2];
    int   smooth_iters = (argc > 3) ? std::atoi(argv[3]) : 5;
    float target_ratio = (argc > 4) ? std::atof(argv[4]) : 0.5f;

    std::fprintf(stderr, "=== Mesh Optimizer ===\n");
    std::fprintf(stderr, "Input:          %s\n", infile);
    std::fprintf(stderr, "Output:         %s\n", outfile);
    std::fprintf(stderr, "Smooth iters:   %d\n", smooth_iters);
    std::fprintf(stderr, "Target ratio:   %.2f\n", target_ratio);
#ifdef _OPENMP
    std::fprintf(stderr, "OpenMP threads: %d\n", omp_get_max_threads());
#endif
    std::fprintf(stderr, "\n");

    // 1. Read
    std::fprintf(stderr, "[1/3] Reading STL ...\n");
    Mesh m = read_stl(infile);

    // 2. Smooth
    if (smooth_iters > 0) {
        std::fprintf(stderr, "[2/3] Laplacian smoothing (%d passes) ...\n", smooth_iters);
        Timer t("smooth total");
        laplacian_smooth(m, smooth_iters);
    } else {
        std::fprintf(stderr, "[2/3] Smoothing skipped\n");
    }

    // 3. Decimate
    if (target_ratio < 1.0f) {
        std::fprintf(stderr, "[3/3] Decimating ...\n");
        decimate(m, target_ratio);
    } else {
        std::fprintf(stderr, "[3/3] Decimation skipped\n");
    }

    // 4. Write — detect output format by extension
    std::fprintf(stderr, "\nWriting output ...\n");
    if (has_extension(outfile, "h5mesh") || has_extension(outfile, "h5")) {
#if HAVE_HDF5
        write_h5mesh(outfile, m);
#else
        std::fprintf(stderr,
            "Error: HDF5 output requested but not compiled with HDF5 support.\n"
            "Rebuild with: g++ -O2 -std=c++17 -fopenmp "
            "-I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5 "
            "-o mesh_optimize mesh_optimize.cpp\n");
        return 1;
#endif
    } else {
        write_stl(outfile, m);
    }
    std::fprintf(stderr, "Done.\n");

    return 0;
}
