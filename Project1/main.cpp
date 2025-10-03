#include <glad/glad.h>
#include <glad/glad.c> // NOWE — musi byæ PRZED GLFW
#include <GLFW/glfw3.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <thread>           // std::thread, hardware_concurrency
#include <algorithm>        // std::sort, std::unique, std::min/std::max
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#define M_PI 3.1415

#include <cstdlib> // wymagane dla exit()

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void renderScene();
//void drawCuboid(float width, float height, float depth);
void drawCuboidTransparentSorted();
int printOversizedTriangles(float maxArea);

struct node {
    glm::vec3 position = glm::vec3(0, 0, 0); // Domyslnie ustawione na (0,0,0)
    glm::vec3 velocity = glm::vec3(0, 0, 0);
    float energy = 0;
};

std::vector<node> nodes;
std::vector<node> mic_nodes;

struct Cuboid_dimensions {
    float width = 20.0f;
    float height = 12.0f;
    float depth = 20.0f;
};
Cuboid_dimensions Cube;

//MIKROFON
struct Microphone_pos{
    float mic_x;
    float mic_y;
    float mic_z;
};
Microphone_pos Mic_pos;


struct Triangle {
    int indices[3];
};
std::vector<Triangle> triangles;
std::vector<Triangle> microphone;

void drawIcosahedron(float radius, std::vector<Triangle>);
void drawMicrophone(float mic_radius, std::vector<Triangle>);
void checkMicrophone(float mic_radius, std::vector<node> wave, std::vector<node> mic);
int pruneSlowNodes(float minSpeed);


// unikalny klucz krawêdzi (mniejszy indeks najpierw)
static inline uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (uint64_t(a) << 32) | uint32_t(b);
}

// rzut na sferê (promieñ r)
static inline glm::vec3 normalize_to_radius(const glm::vec3& p, float r) {
    float len = glm::length(p);
    if (len == 0.0f) return glm::vec3(r, 0, 0);
    return (p / len) * r;
}

// zwraca indeks œrodkowego wierzcho³ka (z cache'owaniem, by unikn¹æ duplikatów)
static int midpoint_index(int i0, int i1,
    std::vector<glm::vec3>& verts,
    std::unordered_map<uint64_t, int>& cache,
    float radius)
{
    uint64_t key = edge_key(i0, i1);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    glm::vec3 m = 0.5f * (verts[i0] + verts[i1]);
    m = normalize_to_radius(m, radius);       // PROJEKCJA NA SFERÊ
    int idx = (int)verts.size();
    verts.push_back(m);
    cache.emplace(key, idx);
    return idx;
}

void updatePhysics(float dt) {
    const float cuboidHalfWidth = Cube.width / 2.0f;
    const float cuboidHalfHeight = Cube.height / 2.0f;
    const float cuboidHalfDepth = Cube.depth / 2.0f;
    const float elasticity = 0.8f; // wspoolczynnik sprezystosci (0.8 = 80% energii zachowane)
    

    for (auto& node : nodes) {
        // Aktualizacja pozycji
        node.position = node.position + node.velocity * dt;
        node.energy = node.energy * 0.999;

        // Sprawdzanie kolizji ze scianami prostopadloscianu i odbicia
        

        if (node.position.y <= -cuboidHalfHeight || node.position.y >= cuboidHalfHeight)
        {
            node.velocity.y = -node.velocity.y * elasticity;
            node.position.y = (node.position.y < 0) ? -cuboidHalfHeight + 0.001f : cuboidHalfHeight - 0.001f;
        }

        if (node.position.x <= -cuboidHalfWidth || node.position.x >= cuboidHalfWidth)
        {
            node.velocity.x = -node.velocity.x * elasticity;
            // Korekta pozycji aby nie utknac w scianie
            node.position.x = (node.position.x < 0) ? -cuboidHalfWidth + 0.001f : cuboidHalfWidth - 0.001f;
        }

        if (node.position.z <= -cuboidHalfDepth || node.position.z >= cuboidHalfDepth)
        {
            node.velocity.z = -node.velocity.z * elasticity;
            node.position.z = (node.position.z < 0) ? -cuboidHalfDepth + 0.001f : cuboidHalfDepth - 0.001f;
        }
    }

}

float calculateTriangleArea(int a, int b, int c) {
    glm::vec3 ab = nodes[b].position - nodes[a].position;
    glm::vec3 ac = nodes[c].position - nodes[a].position;
    glm::vec3 cross = glm::cross(ab, ac);
    return 0.5f * glm::length(cross); //funkcja prawdziwego pola
    //return glm::dot(cross, cross); // pole bez pierwiastka czyli (2*area)^2
}

int addMidpoint(int a, int b) {
    node midpoint;
    midpoint.position = (nodes[a].position + nodes[b].position) * 0.5f;
    midpoint.velocity = (nodes[a].velocity + nodes[b].velocity) * 0.5f;
    //midpoint.velocity = (nodes[a].velocity + nodes[b].velocity) * 0.5f;

    nodes.push_back(midpoint);
    return (int)(nodes.size() - 1);
}

// do podzialu sfery na klatki
// ==== PERSISTENT STATE dla rafinowania „po kawa³ku” ====
static size_t gRefineCursor = 0; // gdzie skoñczyliœmy ostatnio
static std::unordered_map<uint64_t, int> gEdgeMidCache; // edge -> midpoint

// ====== POMOCNICZE ======
struct Edge { int a, b; }; // zawsze a<b
static inline Edge make_edge(int u, int v) { if (u > v) std::swap(u, v); return { u,v }; }
static inline bool edge_less(const Edge& x, const Edge& y) {
    return (x.a < y.a) || (x.a == y.a && x.b < y.b);
}
static inline bool edge_eq(const Edge& x, const Edge& y) {
    return x.a == y.a && x.b == y.b;
}

static int midpoint_index_nodes_cached(int i0, int i1) {
    uint64_t key = edge_key(i0, i1);           // masz ju¿ edge_key(i,j)
    auto it = gEdgeMidCache.find(key);
    if (it != gEdgeMidCache.end()) return it->second;
    int idx = addMidpoint(i0, i1);      // Twój kod midpointu (pozycja/prêdkoœæ)
    gEdgeMidCache.emplace(key, idx);
    return idx;
}

// Zwraca indeks midpointu z cache; jeœli nie istnieje — tworzy przez addMidpoint(...)
static int midpoint_index_nodes(int i0, int i1,
    std::unordered_map<uint64_t, int>& cache)
{
    uint64_t key = edge_key(i0, i1);     // edge_key ju¿ masz w pliku
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;

    int idx = addMidpoint(i0, i1); // œrednia pozycji i prêdkoœci — jak u Ciebie
    cache.emplace(key, idx);
    return idx;
}

void refineIcosahedron(std::vector<node>& nodes,
    float maxArea)
{
    std::vector<Triangle> newTriangles;
    newTriangles.reserve(triangles.size() * 4);

    // cache midpointów na czas ca³ej iteracji podzia³u
    std::unordered_map<uint64_t, int> midCache;
    midCache.reserve(triangles.size() * 3);

    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.width * 0.5f;
    const float safetyMargin = 0.2f;

    // Porównujemy bez sqrt: |cross|^2 > (2*maxArea)^2
    const float area2_threshold = 4.0f * maxArea * maxArea;

    for (const auto& tri : triangles) {
        const int a = tri.indices[0];
        const int b = tri.indices[1];
        const int c = tri.indices[2];

        // Lokalnie trzymaj pozycje (mniej dostêpów do pamiêci)
        const glm::vec3& pa = nodes[a].position;
        const glm::vec3& pb = nodes[b].position;
        const glm::vec3& pc = nodes[c].position;

        // 1) Tanie sprawdzenie nearWall — szybciej od razu odpuœciæ
        const bool nearWall =
            (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

        bool shouldSplit = false;
        if (!nearWall) {
            // 2) Licz „pole” bez sqrt: |cross|^2 vs threshold
            const glm::vec3 ab = pb - pa;
            const glm::vec3 ac = pc - pa;
            const glm::vec3 cr = glm::cross(ab, ac);
            const float     cr2 = glm::dot(cr, cr);
            shouldSplit = (cr2 > area2_threshold);
        }

        if (shouldSplit) {
            // 3) Midpointy z cache — brak duplikatów przy krawêdziach wspó³dzielonych
            const int ab_i = midpoint_index_nodes(a, b, midCache);
            const int bc_i = midpoint_index_nodes(b, c, midCache);
            const int ca_i = midpoint_index_nodes(c, a, midCache);

            newTriangles.push_back({ {a,  ab_i, ca_i} });
            newTriangles.push_back({ {b,  bc_i, ab_i} });
            newTriangles.push_back({ {c,  ca_i, bc_i} });
            newTriangles.push_back({ {ab_i, bc_i, ca_i} });
        }
        else {
            newTriangles.push_back(tri);
        }
    }

    triangles.swap(newTriangles);
}

// ==== NOWA WERSJA: rafinowanie bud¿etowe ====
bool refineIcosahedron_chunked(float maxArea,
    size_t triBudget /* ile trójk¹tów obrabiamy na wywo³anie */)
{
    if (triangles.empty() || triBudget == 0) return false;

    // Jeœli od ostatniego razu siatka zosta³a zresetowana – wyzeruj stan
    if (gRefineCursor > triangles.size()) {
        gRefineCursor = 0;
        gEdgeMidCache.clear();
    }

    // Sta³e i próg bez sqrt: |cross|^2 > (2*maxArea)^2
    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.depth * 0.5f;
    const float safetyMargin = 0.2f;
    const float area2_threshold = 4.0f * maxArea * maxArea;

    // Pracujemy tylko na partii trójk¹tów, które ISTNIA£Y na wejœciu
    const size_t size0 = triangles.size();
    size_t start = (gRefineCursor < size0 ? gRefineCursor : 0);
    size_t end = std::min(start + triBudget, size0);

    bool changed = false;

    for (size_t i = start; i < end; ++i) {
        // UWAGA: nie trzymaj referencji do triangles[i] — bêdziemy push_back’owaæ.
        //Triangle t = triangles[i];
        int a = triangles[i].indices[0], b = triangles[i].indices[1], c = triangles[i].indices[2]; // tutaj zamiast triangles[i] bylo t.

        const glm::vec3& pa = nodes[a].position;
        const glm::vec3& pb = nodes[b].position;
        const glm::vec3& pc = nodes[c].position;

        // Tanie sprawdzenie „nearWall” — jeœli blisko œciany, nie dzielimy
        const bool nearWall =
            (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

        if (!nearWall) {
            // Pole bez sqrt: |(pb-pa) x (pc-pa)|^2
            const glm::vec3 ab = pb - pa;
            const glm::vec3 ac = pc - pa;
            const glm::vec3 cr = glm::cross(ab, ac);
            const float     cr2 = glm::dot(cr, cr);

            if (cr2 > area2_threshold) {
                // Zast¹p trójk¹t i do³ó¿ 3 nowe — in-place, bez kopiowania ca³ej tablicy
                const int ab_i = midpoint_index_nodes_cached(a, b);
                const int bc_i = midpoint_index_nodes_cached(b, c);
                const int ca_i = midpoint_index_nodes_cached(c, a);

                triangles[i] = { { a,  ab_i, ca_i } };              // zamiana bie¿¹cego
                triangles.push_back({ { b,  bc_i, ab_i } });        // 3 nowe
                triangles.push_back({ { c,  ca_i, bc_i } });
                triangles.push_back({ { ab_i, bc_i, ca_i } });

                changed = true;
            }
        }
    }

    // Przesuñ kursor; po dojœciu do koñca partii – owiñ na pocz¹tek
    gRefineCursor = end;
    if (gRefineCursor >= size0) gRefineCursor = 0;

    return changed;
}

// ====== MULTI-THREAD CHUNKED REFINEMENT ======
bool refineIcosahedron_chunked_mt(float maxArea,
    size_t triBudget,          // ile trójk¹tów obrabiamy w tej klatce
    int   threadCount = std::thread::hardware_concurrency())
{
    if (triangles.empty() || triBudget == 0) return false;
    if (threadCount <= 0) threadCount = 1;

    // jeœli siatka by³a resetowana
    if (gRefineCursor > triangles.size()) {
        gRefineCursor = 0;
        gEdgeMidCache.clear();
    }

    // próg bez sqrt: |cross|^2 > (2*maxArea)^2
    const float area2_threshold = 4.0f * maxArea * maxArea;

    // sta³e kolizyjne (jak u Ciebie)
    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.depth * 0.5f;
    const float safetyMargin = 0.2f;

    // pracujemy wy³¹cznie na „starym” prefiksie tablicy trójk¹tów
    const size_t size0 = triangles.size();
    size_t start = (gRefineCursor < size0 ? gRefineCursor : 0);
    size_t end = std::min(start + triBudget, size0);
    const size_t N = (end > start ? end - start : 0);
    if (N == 0) return false;

    threadCount = std::min<int>(threadCount, (int)N);

    // --- FAZA 1: RÓWNOLEGLE zbierz „do podzia³u” + krawêdzie (TLS) ---
    std::vector<std::vector<int>>   tl_split(threadCount);
    std::vector<std::vector<Edge>>  tl_edges(threadCount);

    auto worker_scan = [&](int tid) {
        size_t chunk = (N + threadCount - 1) / threadCount;
        size_t i0 = start + tid * chunk;
        size_t i1 = std::min(end, i0 + chunk);
        auto& split = tl_split[tid];
        auto& edges = tl_edges[tid];
        split.reserve((i1 > i0) ? (i1 - i0) / 2 : 0);
        edges.reserve((i1 > i0) ? 3 * (i1 - i0) / 2 : 0);

        for (size_t i = i0; i < i1; ++i) {
            Triangle t = triangles[i]; // kopiuj (nie referencja; póŸniej bêdziemy pushowaæ do triangles)
            int a = t.indices[0], b = t.indices[1], c = t.indices[2];

            const glm::vec3& pa = nodes[a].position;
            const glm::vec3& pb = nodes[b].position;
            const glm::vec3& pc = nodes[c].position;

            const bool nearWall =
                (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
                (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
                (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
                (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
                (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

            if (nearWall) continue;

            const glm::vec3 ab = pb - pa;
            const glm::vec3 ac = pc - pa;
            const glm::vec3 cr = glm::cross(ab, ac);
            const float     cr2 = glm::dot(cr, cr);

            if (cr2 > area2_threshold) {
                split.push_back((int)i);
                edges.push_back(make_edge(a, b));
                edges.push_back(make_edge(b, c));
                edges.push_back(make_edge(c, a));
            }
        }
        };

    std::vector<std::thread> threads;
    threads.reserve(threadCount);
    for (int t = 0; t < threadCount; ++t) threads.emplace_back(worker_scan, t);
    for (auto& th : threads) th.join();
    threads.clear();

    // scalenie wyników fazy 1
    std::vector<int> splitIdx;
    std::vector<Edge> edges;
    {
        size_t totalSplit = 0, totalEdges = 0;
        for (int t = 0; t < threadCount; ++t) { totalSplit += tl_split[t].size(); totalEdges += tl_edges[t].size(); }
        if (totalSplit == 0) { // nic do roboty
            gRefineCursor = end < size0 ? end : 0;
            return false;
        }
        splitIdx.reserve(totalSplit);
        edges.reserve(totalEdges);
        for (int t = 0; t < threadCount; ++t) {
            splitIdx.insert(splitIdx.end(), tl_split[t].begin(), tl_split[t].end());
            edges.insert(edges.end(), tl_edges[t].begin(), tl_edges[t].end());
            tl_split[t].clear(); tl_edges[t].clear();
        }
    }

    // unikalizacja krawêdzi
    std::sort(edges.begin(), edges.end(), edge_less);
    edges.erase(std::unique(edges.begin(), edges.end(), edge_eq), edges.end());

    // --- FAZA 2: SEKWENCYJNIE utwórz midpointy (u¿yj globalnego cache) ---
    std::vector<int> edgeMidIdx(edges.size());
    nodes.reserve(nodes.size() + edges.size());          // minimalizacja realokacji
    for (size_t i = 0; i < edges.size(); ++i) {
        edgeMidIdx[i] = midpoint_index_nodes_cached(edges[i].a, edges[i].b);
    }
    auto edge_lookup = [&](int u, int v)->int {
        if (u > v) std::swap(u, v);
        Edge key{ u,v };
        auto it = std::lower_bound(edges.begin(), edges.end(), key, edge_less);
        // przy poprawnym zbiorze powinien istnieæ
        return edgeMidIdx[size_t(it - edges.begin())];
        };

    // --- FAZA 3: RÓWNOLEGLE zbuduj nowe trójk¹ty do buforów per-w¹tek ---
    std::vector<std::vector<std::pair<int, Triangle>>> tl_replace(threadCount); // (index, tri)
    std::vector<std::vector<Triangle>>                tl_append(threadCount);  // 3 szt./split

    auto worker_build = [&](int tid) {
        size_t M = splitIdx.size();
        size_t chunk = (M + threadCount - 1) / threadCount;
        size_t i0 = tid * chunk;
        size_t i1 = std::min(M, i0 + chunk);
        auto& rep = tl_replace[tid];
        auto& app = tl_append[tid];
        rep.reserve((i1 > i0) ? (i1 - i0) : 0);
        app.reserve((i1 > i0) ? 3 * (i1 - i0) : 0);

        for (size_t k = i0; k < i1; ++k) {
            int ti = splitIdx[k];
            Triangle t = triangles[ti];
            int a = t.indices[0], b = t.indices[1], c = t.indices[2];

            int ab = edge_lookup(a, b);
            int bc = edge_lookup(b, c);
            int ca = edge_lookup(c, a);

            // zamieniamy bie¿¹cy
            rep.emplace_back(ti, Triangle{ {a,  ab, ca} });
            // 3 nowe
            app.push_back({ {b,  bc, ab} });
            app.push_back({ {c,  ca, bc} });
            app.push_back({ {ab, bc, ca} });
        }
        };

    for (int t = 0; t < threadCount; ++t) threads.emplace_back(worker_build, t);
    for (auto& th : threads) th.join();

    // --- FAZA 4: JEDNYM KROKIEM zapisz wynik do 'triangles' ---
    // (po fazie 3 nie czytamy ju¿ 'triangles', wiêc mo¿emy rezerwowaæ i pushowaæ)
    size_t toAppend = 0;
    for (int t = 0; t < threadCount; ++t) toAppend += tl_append[t].size();
    triangles.reserve(triangles.size() + toAppend);

    // podmiany
    for (int t = 0; t < threadCount; ++t) {
        for (auto& p : tl_replace[t]) {
            triangles[p.first] = p.second;
        }
    }
    // dopisywanie 3 nowych
    for (int t = 0; t < threadCount; ++t) {
        triangles.insert(triangles.end(), tl_append[t].begin(), tl_append[t].end());
    }

    // przesuwamy kursor (zawsze wzglêdem dawnego size0)
    gRefineCursor = end;
    if (gRefineCursor >= size0) gRefineCursor = 0;

    return true;
}


const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

bool firstMouse = true;
bool mousePressed = false;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
float yaw = -90.0f;
float pitch = 0.0f;
float fov = 45.0f;
float deltaTime = 0.0f;
float lastFrame = 0.0f;
bool first = true;

//do przyspieszenia VBO
static GLuint gVBO = 0;
static GLuint gIBO = 0;
static GLsizei gIndexCount = 0;
static bool gMeshDirty = true;           // trzeba odbudowaæ bufory (zmiana triangles/nodes count)
static size_t gLastV = 0, gLastI = 0;    // do detekcji zmian topologii

//proba
// VBO/IBO/VAO
static GLuint gVAO = 0;
static GLuint gVboPos = 0;      // tylko pozycje (glm::vec3)
static GLuint gIbo = 0;

static GLsizeiptr gPosBytes = 0;
static GLsizeiptr gIboBytes = 0;

// persistent mapping (opcjonalnie)
static bool       gPosPersistent = false;
static glm::vec3* gPosPtr = nullptr; // wskazanie na zmapowany bufor pozycji

// do FPS
int frameCount = 0;
double previousTime = 0.0;
double fps = 0.0;

void calculateFPS() {
    double currentTime = glfwGetTime();
    frameCount++;

    if (currentTime - previousTime >= 1.0) {
        fps = frameCount / (currentTime - previousTime);
        std::cout << "FPS: " << fps << std::endl;

        frameCount = 0;
        previousTime = currentTime;
        printOversizedTriangles(0.05f);
    }
}

float radius = 0.2f;
float mic_radius = 0.1f;
const float H_ANGLE = M_PI / 180 * 72; // 72 stopni w radianach
const float V_ANGLE = atanf(1.0f / 2); // Kat wierzcholka

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            mousePressed = true;
        }
        else if (action == GLFW_RELEASE) {
            mousePressed = false;
            firstMouse = true;
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    // Zmiana odleglosci kamery za pomoca scrolla
    cameraPos += cameraFront * static_cast<float>(yoffset) * 0.5f;
}

static void glfwErrorCallback(int code, const char* desc) { fprintf(stderr, "GLFW[%d]: %s\n", code, desc); }

int main()
{
    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Camera Scene", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // to pod tym dodane do VBO 
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "GLAD fail\n"; return -1;
    }
    std::cout << "GL " << GLVersion.major << "." << GLVersion.minor << "\n";
    //glEnable(GL_DEPTH_TEST);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    //--POCZATKOWA POZYCJA MIKROFONU---
    Mic_pos.mic_x = 3.5f;
    Mic_pos.mic_y = 2.0f;
    Mic_pos.mic_z = 1.0f;

    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);

        glClearColor(0.5f, 0.5f, 0.5f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(fov), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(&projection[0][0]);

        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(&view[0][0]);

        renderScene();
        //ZMIEN POZYCJE MIKROFONU
        Mic_pos.mic_x = Mic_pos.mic_x - 1.0f;

        // Oblicz FPS
        calculateFPS();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void processInput(GLFWwindow* window)
{
    float cameraSpeed = 2.5f * deltaTime;

    // Zmiana wysokosci kamery za pomoca klawiszy W i S
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos.y += cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos.y -= cameraSpeed;

    // Poruszanie kamera w plaszczyznie XZ za pomoca klawiszy A i D
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (!mousePressed)
        return;

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    if (pitch > 89.0f)
        pitch = 89.0f;
    if (pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(front);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void buildSphereBuffers(bool dynamic = true) {
    if (!GLAD_GL_VERSION_1_5) { std::cerr << "VBO niedostêpne (GL < 1.5)\n"; return; }

    // --- VBO (TERAZ: ca³y array 'nodes', bez przepakowywania do vec<vec3>) ---
    if (!gVBO) glGenBuffers(1, &gVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gVBO);

    const GLsizeiptr vbSize = (GLsizeiptr)(nodes.size() * sizeof(node));
    const GLenum usage = dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW;

    glBufferData(GL_ARRAY_BUFFER, vbSize, nodes.empty() ? nullptr : (const void*)nodes.data(), usage);

    // Atrybut pozycji: 3 floaty, stride = sizeof(node), offset = offsetof(node, position)
    glEnableClientState(GL_VERTEX_ARRAY); // compatibility profile
    glVertexPointer(3, GL_FLOAT, (GLsizei)sizeof(node), (const void*)offsetof(node, position));

    // --- IBO (bez tworzenia wektora 'indices') ---
    if (!gIBO) glGenBuffers(1, &gIBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIBO);

    const GLsizeiptr ibSize = (GLsizeiptr)(triangles.size() * 3 * sizeof(unsigned int));
    const void* idxSrc = triangles.empty() ? nullptr : (const void*)&triangles[0].indices[0];

    glBufferData(GL_ELEMENT_ARRAY_BUFFER, ibSize, idxSrc, GL_STATIC_DRAW);

    gIndexCount = (GLsizei)(triangles.size() * 3);

    // porz¹dek
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
}


void updateSpherePositions() {
    if (!gVBO) return;

    // Uwaga: wysy³amy CA£¥ strukturê node (nie tylko pozycjê)
    // — ale to eliminuje pêtlê kopiuj¹c¹ pozycje po CPU.
    glBindBuffer(GL_ARRAY_BUFFER, gVBO);
    const GLsizeiptr vbSize = (GLsizeiptr)(nodes.size() * sizeof(node));

    // orphan + jednorazowy transfer
    glBufferData(GL_ARRAY_BUFFER, vbSize, nullptr, GL_DYNAMIC_DRAW);
    if (!nodes.empty())
        glBufferSubData(GL_ARRAY_BUFFER, 0, vbSize, (const void*)nodes.data());

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void drawSphereWithBuffers()
{
    if (!gVBO || !gIBO || gIndexCount == 0) return;

    glBindBuffer(GL_ARRAY_BUFFER, gVBO);
    glEnableClientState(GL_VERTEX_ARRAY);                           // compat profile
    glVertexPointer(3, GL_FLOAT, (GLsizei)sizeof(node),
        (const void*)offsetof(node, position));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIBO);

    // 1) Wype³nienie – z polygon offset, ¿eby obrys nie „miga³”
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0f, 1.0f);          // odsuñ fill w g³¹b z-bufferu

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
    glDrawElements(GL_TRIANGLES, gIndexCount, GL_UNSIGNED_INT, (void*)0);

    glDisable(GL_POLYGON_OFFSET_FILL);

    // 2) Obrys (wireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(1.0f);                    // mo¿esz podbiæ np. do 2.0
    glColor4f(0.0f, 0.0f, 0.0f, 0.6f);
    glDrawElements(GL_TRIANGLES, gIndexCount, GL_UNSIGNED_INT, (void*)0);

    // 3) Sprz¹tanie
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void renderScene()
{
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    float dt = 0.01;
    //static std::vector<Triangle> triangles;
    if (first) {
        nodes.clear();
        triangles.clear();

        //gorny wierzcholek
        node buf;
        buf.position = glm::vec3(Mic_pos.mic_x, Mic_pos.mic_y + mic_radius, Mic_pos.mic_z);
        mic_nodes.push_back(buf);

        //boczne wierzcholki
        for (int i = 1; i <= 5; ++i) {
            buf.position = glm::vec3(Mic_pos.mic_x + mic_radius * cos(V_ANGLE) * cos(i * H_ANGLE), Mic_pos.mic_y + mic_radius * sin(V_ANGLE), Mic_pos.mic_z + mic_radius * cos(V_ANGLE) * sin(i * H_ANGLE));
            mic_nodes.push_back(buf);
        }
        for (int i = 6; i <= 10; ++i) {
            buf.position = glm::vec3(Mic_pos.mic_x + mic_radius * cos(V_ANGLE) * cos((i + 0.5f) * H_ANGLE), Mic_pos.mic_y + -mic_radius * sin(V_ANGLE), Mic_pos.mic_z + mic_radius * cos(V_ANGLE) * sin((i + 0.5f) * H_ANGLE));
            mic_nodes.push_back(buf);
        }

        //dolny wierzcholek
        buf.position = glm::vec3(Mic_pos.mic_x, Mic_pos.mic_y - mic_radius, Mic_pos.mic_z);
        mic_nodes.push_back(buf);

        const int initTris[20][3] = {
            {0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, {0,5,1},
            {11,6,7}, {11,7,8}, {11,8,9}, {11,9,10}, {11,10,6},
            {1,2,6}, {2,3,7}, {3,4,8}, {4,5,9}, {5,1,10},
            {6,7,2}, {7,8,3}, {8,9,4}, {9,10,5}, {10,6,1}
        };

        for (int i = 0; i < 20; ++i) {
            microphone.push_back({ {initTris[i][0], initTris[i][1], initTris[i][2]} });
        }

        // --- Parametr gêstoœci: ka¿dy poziom ×4 liczba trójk¹tów ---
        constexpr int SUBDIV = 2; // 2 optymalnie, wiecej laguje

        // --- 12 wierzcho³ków  na sferze o promieniu 'radius' ---
        const float t = (1.0f + std::sqrt(5.0f)) * 0.5f; // z³ota proporcja
        std::vector<glm::vec3> verts = {
            glm::normalize(glm::vec3(-1,  t,  0)),
            glm::normalize(glm::vec3(1,  t,  0)),
            glm::normalize(glm::vec3(-1, -t,  0)),
            glm::normalize(glm::vec3(1, -t,  0)),

            glm::normalize(glm::vec3(0, -1,  t)),
            glm::normalize(glm::vec3(0,  1,  t)),
            glm::normalize(glm::vec3(0, -1, -t)),
            glm::normalize(glm::vec3(0,  1, -t)),

            glm::normalize(glm::vec3(t,  0, -1)),
            glm::normalize(glm::vec3(t,  0,  1)),
            glm::normalize(glm::vec3(-t,  0, -1)),
            glm::normalize(glm::vec3(-t,  0,  1)),
        };
        for (auto& v : verts) v *= radius; // na promieñ 'radius'

        // --- 20 œcian (CCW patrz¹c z zewn¹trz) ---
        std::vector<glm::ivec3> faces = {
            {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
            {1,5,9},  {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
            {3,9,4},  {3,4,2},  {3,2,6},  {3,6,8},  {3,8,9},
            {4,9,5},  {2,4,11}, {6,2,10}, {8,6,7},  {9,8,1}
        };

        // --- Subdivision ---
        for (int s = 0; s < SUBDIV; ++s) {
            std::unordered_map<uint64_t, int> cache;
            std::vector<glm::ivec3> new_faces;
            new_faces.reserve(faces.size() * 4);

            for (const auto& f : faces) {
                int i0 = f.x, i1 = f.y, i2 = f.z;

                int a = midpoint_index(i0, i1, verts, cache, radius);
                int b = midpoint_index(i1, i2, verts, cache, radius);
                int c = midpoint_index(i2, i0, verts, cache, radius);

                // 4 nowe trójk¹ty (CCW)
                new_faces.push_back({ i0, a, c });
                new_faces.push_back({ i1, b, a });
                new_faces.push_back({ i2, c, b });
                new_faces.push_back({ a,  b, c });
            }
            faces.swap(new_faces);
        }

        nodes.reserve(verts.size());
        for (const auto& p : verts) {
            node nd;
            nd.position = p;
            nd.velocity = nd.position * 1.0f; 
            nodes.push_back(nd);
        }

        triangles.reserve(faces.size());
        for (const auto& f : faces) {
            triangles.push_back({ { f.x, f.y, f.z } });
            //break;
        }

        first = false;
        buildSphereBuffers(/*dynamic=*/true);
    }

    //MIKROFON TUTAJ
    drawMicrophone(mic_radius, microphone);
    checkMicrophone(mic_radius, nodes, mic_nodes);

    //buildSphereBuffers(nodes, triangles, /*dynamic=*/true); //bledne bo caly czas sie wykonywalo

    //std::cout << "Ilosc punktow:" << nodes.size() << std::endl;
    //if (nodes.size() >= 10000) exit(0);
    //std::cout << "Ilosc scian:" << triangles.size() << std::endl;
    /*for (int i = 0; i < nodes.size(); i++)
    {
        std::cout << "X:" << nodes[i].position.x << std::endl;
        std::cout << "Y:" << nodes[i].position.y << std::endl;
        std::cout << "Z:" << nodes[i].position.z << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "koniec" << std::endl;*/

    //static int frameCount = 0;
    //if (frameCount % 10 == 0) {
    //    refineIcosahedron(nodes, 0.05f);
    //    gMeshDirty = true;   // <-- DODAJ
    //}
    //frameCount++;

    static int frameCount = 0;
    if ((frameCount % 2) == 0) {
        size_t budget = std::min<size_t>(4000, std::max<size_t>(1, triangles.size() / 10));
        budget = triangles.size() / 12;
        int threads = std::max(1u, std::thread::hardware_concurrency());

        if (refineIcosahedron_chunked_mt(0.05f, budget,threads)) {
            gMeshDirty = true; // przebuduj bufory tylko gdy zasz³a zmiana
        }

        /*refineIcosahedron(nodes, 0.05f);
        gMeshDirty = true;
        */
    }
    frameCount++;

    //updatePhysics(dt);
    /*Cuboid_dimensions Cube;
    Cube.width = 80.0f;
    Cube.height = 60.0f;
    Cube.depth = 100.0f;*/
    updatePhysics(dt);
    //int removed = pruneSlowNodes(nodes, /*minSpeed=*/2.00f); // m/s (dobierz)

    // 1) odbuduj bufory tylko gdy trzeba
    if (gMeshDirty || nodes.size() != gLastV || triangles.size() != gLastI) {
        buildSphereBuffers(/*dynamic=*/true);
        gLastV = nodes.size();
        gLastI = triangles.size();
        gMeshDirty = false;
    }
    else {
        // 2) w przeciwnym razie tylko podmieñ pozycje
        updateSpherePositions();
    }

    // 3) rysuj TYLKO z VBO/IBO
    drawSphereWithBuffers();

    
    //kulka
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    //drawIcosahedron(0.2f, nodes, triangles);
    //drawSphereWithBuffers();
    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Basen
    glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
    //drawCuboid(Cube.width, Cube.height, Cube.depth);
    drawCuboidTransparentSorted();
}

int printOversizedTriangles(float maxArea) {
    // te same sta³e co w refine
    const float cuboidHalfWidth = Cube.width * 0.5f;
    const float cuboidHalfHeight = Cube.height * 0.5f;
    const float cuboidHalfDepth = Cube.depth * 0.5f;
    const float safetyMargin = 0.2f;

    // porównujemy |cross|^2 > (2*maxArea)^2
    const float thr2 = 4.0f * maxArea * maxArea;

    int count = 0;
    for (const auto& t : triangles) {
        int a = t.indices[0], b = t.indices[1], c = t.indices[2];
        const glm::vec3& pa = nodes[a].position;
        const glm::vec3& pb = nodes[b].position;
        const glm::vec3& pc = nodes[c].position;

        // pomijamy trójk¹ty blisko œciany (jak w refine)
        const bool nearWall =
            (std::fabs(pa.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pa.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pa.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pb.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pb.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pb.z) > cuboidHalfDepth - safetyMargin) ||
            (std::fabs(pc.x) > cuboidHalfWidth - safetyMargin) ||
            (std::fabs(pc.y) > cuboidHalfHeight - safetyMargin) ||
            (std::fabs(pc.z) > cuboidHalfDepth - safetyMargin);

        if (nearWall) continue;

        // pole bez sqrt: |(pb-pa) x (pc-pa)|^2
        const glm::vec3 cr = glm::cross(pb - pa, pc - pa);
        const float cr2 = glm::dot(cr, cr);

        if (cr2 > thr2) ++count;
    }
    float procent = (float)count / (float)triangles.size() * 100.0f;
    std::cout << "Oversized (splittable) triangles: "
        << triangles.size() << "\n";
    return count;
}

// Usuwa wêz³y o |velocity| < minSpeed i remapuje trójk¹ty.
// Zwraca ile wêz³ów usuniêto.
int pruneSlowNodes(float minSpeed) {
    const float thr2 = minSpeed * minSpeed;
    const size_t N = nodes.size();

    // mapowanie stary_index -> nowy_index (=-1 gdy usuniêty)
    std::vector<int> remap(N, -1);
    std::vector<node> kept;
    kept.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        const glm::vec3& v = nodes[i].velocity;
        const float v2 = glm::dot(v, v);
        if (v2 >= thr2) {
            remap[i] = (int)kept.size();
            kept.push_back(nodes[i]);
        }
    }

    if (kept.size() == N) {
        return 0; // nic nie usuniêto
    }

    // przebuduj trójk¹ty — zostaw tylko te, które w ca³oœci przetrwa³y
    std::vector<Triangle> newTris;
    newTris.reserve(triangles.size());
    for (const auto& t : triangles) {
        int a = t.indices[0], b = t.indices[1], c = t.indices[2];
        int na = remap[a], nb = remap[b], nc = remap[c];
        if (na >= 0 && nb >= 0 && nc >= 0) {
            newTris.push_back({ { na, nb, nc } });
        }
        // inaczej: trójk¹t wypad³, bo dotyka³ usuniêtego wêz³a
    }

    nodes.swap(kept);
    triangles.swap(newTris);

    // uniewa¿nij stan rafinowania „chunked”
    gRefineCursor = 0;
    gEdgeMidCache.clear();

    // poproœ o odbudowê VBO
    gMeshDirty = true;

    return (int)(N - nodes.size());
}

void drawCuboidTransparentSorted() {
    float halfWidth = Cube.width / 2.0f;
    float halfHeight = Cube.height / 2.0f;
    float halfDepth = Cube.depth / 2.0f;

    glm::vec3 vertices[] = {
        { -halfWidth, -halfHeight,  halfDepth },
        {  halfWidth, -halfHeight,  halfDepth },
        {  halfWidth,  halfHeight,  halfDepth },
        { -halfWidth,  halfHeight,  halfDepth },
        { -halfWidth, -halfHeight, -halfDepth },
        {  halfWidth, -halfHeight, -halfDepth },
        {  halfWidth,  halfHeight, -halfDepth },
        { -halfWidth,  halfHeight, -halfDepth }
    };

    struct Face {
        int indices[4];
        float distance;
    };

    Face faces[6] = {
        {{0, 1, 2, 3}, 0}, // przód
        {{4, 5, 6, 7}, 0}, // ty³
        {{0, 3, 7, 4}, 0}, // lewo
        {{1, 2, 6, 5}, 0}, // prawo
        {{0, 1, 5, 4}, 0}, // dó³
        {{2, 3, 7, 6}, 0}  // góra
    };

    // Oblicz odleg³oœci œcian od kamery
    for (int i = 0; i < 6; ++i) {
        glm::vec3 center(0.0f);
        for (int j = 0; j < 4; ++j) {
            center += vertices[faces[i].indices[j]];
        }
        center /= 4.0f;
        faces[i].distance = glm::length(center - cameraPos);
    }

    // Sortuj œciany od najdalszej do najbli¿szej
    for (int i = 0; i < 5; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            if (faces[i].distance < faces[j].distance) {
                std::swap(faces[i], faces[j]);
            }
        }
    }

    // Rysuj œciany w odpowiedniej kolejnoœci
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            glVertex3f(vertices[faces[i].indices[j]].x,
                vertices[faces[i].indices[j]].y,
                vertices[faces[i].indices[j]].z);
        }
    }
    glEnd();
}


void drawIcosahedron(float radius, std::vector<Triangle> triangles) {
    // Rysowanie scian (pomarañczowe)
    glColor4f(1.0f, 0.5f, 0.0f, 0.5f); //te 0.5f na koncu nie ma znaczenia jak blend wylaczony
    glBegin(GL_TRIANGLES);
    for (const auto& tri : triangles) {
        for (int j = 0; j < 3; ++j) {
            glVertex3f(nodes[tri.indices[j]].position.x,
                nodes[tri.indices[j]].position.y,
                nodes[tri.indices[j]].position.z);
        }
    }
    glEnd();

    // Rysowanie krawedzi (czarne)
    glColor4f(0.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_LINES);
    for (const auto& tri : triangles) {
        for (int j = 0; j < 3; ++j) {
            int current = tri.indices[j];
            int next = tri.indices[(j + 1) % 3];
            glVertex3f(nodes[current].position.x, nodes[current].position.y, nodes[current].position.z);
            glVertex3f(nodes[next].position.x, nodes[next].position.y, nodes[next].position.z);
        }
    }
    glEnd();
}

void drawMicrophone(float mic_radius, std::vector<Triangle>)
{
    // Rysowanie scian 
    glColor4f(6.0f, 0.4f, 0.8f, 0.75f);
    glBegin(GL_TRIANGLES);
    for (const auto& mic : microphone) {
        for (int j = 0; j < 3; ++j) {
            glVertex3f(mic_nodes[mic.indices[j]].position.x,
                mic_nodes[mic.indices[j]].position.y,
                mic_nodes[mic.indices[j]].position.z);
        }
    }
    glEnd();

    // Rysowanie krawedzi (czarne)
    glColor4f(0.0f, 0.0f, 0.0f, 0.5f);
    glBegin(GL_LINES);
    for (const auto& mic : microphone) {
        for (int j = 0; j < 3; ++j) {
            int current = mic.indices[j];
            int next = mic.indices[(j + 1) % 3];
            glVertex3f(mic_nodes[current].position.x, mic_nodes[current].position.y, mic_nodes[current].position.z);
            glVertex3f(mic_nodes[next].position.x, mic_nodes[next].position.y, mic_nodes[next].position.z);
        }
    }
    glEnd();
}

void checkMicrophone(float mic_radius, std::vector<node> wave, std::vector<node> mic)
{
    for(const auto& wav : wave)
    {
        for (const auto& m : mic)
        {
            /*
            if (abs(wav.position.x - m.position.x) < mic_radius && abs(wav.position.y - m.position.y) < mic_radius && abs(wav.position.z - m.position.z) < mic_radius)
            {
                std::cout << "ODCZYT MIKROFONU" << std::endl;
            }
            */
            if (sqrt( pow((wav.position.x - m.position.x),2) + pow((wav.position.y - m.position.y), 2) + pow((wav.position.z - m.position.z), 2)) < mic_radius)
            {
                std::cout << "ODCZYT MIKROFONU" << std::endl;
            }
        }
    }
}